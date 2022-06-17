"""
    imsrg_main(binfo,Chan1b,Chan2bD,HFobj,dictMono,d9j,HOBs,dict6j,valencespace,Operators,to; core_generator_type="atan",valence_generator_type="shell-model-atan",denominatorDelta=0.0)
    
# Arguments
- `binfo::basedat` struct basedat(nuc::nuclei,sntf::String,hw::Float,emax::Int)
- `Chan1b::chan1b` struct for one-body stuffs
- `Chan2bD::chan2bD` struct for two-body stuffs (e.g., dict to get idx from JPT)
- `HFobj::HamiltonianNormalOrdered` struct HNO, which includes info. of HF solution (HF energy, occupation, f,Gamma,...)
- `dictMono::Dict` dictionary to get Vmonopole
- `d9j` preallocated dictionaries for wigner 9j symbols, which are needed to calculate Rp2
- `HOBs` struct HarmonicOscillatorBrackets 
- `dict6j` preallocated dictionaries for wigner6j symbols (needed in e.g., Pandya transformation)
- `valencespace` to specify valence space  
- `Operators::Vector{String}` non-Hamiltonian operators
- `to` TimerOutput object to measure runtime&memory allocations

# Optional Arguments
- `core_generator_type` only the "atan" is implemented
- `valence_generator_type` only the "shell-model-atan" is implemented
- `denominatorDelta::Float` denominator Delta, which is needed for multi-major shell decoupling
"""
function imsrg_main(binfo,Chan1b,Chan2bD,HFobj,dictsnt,d9j,HOBs,dict6j,valencespace,Operators,to;
                    core_generator_type="atan",valence_generator_type="shell-model-atan")
    dictMono = deepcopy(dictsnt.dictMonopole)
    vsIMSRG = ifelse(valencespace!=[],true,false)
    update_core_in_sps!(binfo,HFobj)
    valencesps = check_valence_space(HFobj,valencespace)
    update_vsspace_chs!(HFobj,valencesps,Chan2bD.Chan2b)

    Chan2b = Chan2bD.Chan2b
    init_dictMonopole!(dictMono,Chan2b)
    IMSRGobj = init_IMSRGobject(HFobj)    
    PandyaObj = prep_PandyaLookup(binfo,HFobj,Chan1b,Chan2bD)

    ##IMSRG(2) calculation for target nucleus
    IMSRGflow(binfo,HFobj,IMSRGobj,PandyaObj,Chan1b,Chan2bD,dictMono,dict6j,core_generator_type,valence_generator_type,to)
    if vsIMSRG
        IMSRGflow(binfo,HFobj,IMSRGobj,PandyaObj,Chan1b,Chan2bD,dictMono,dict6j,core_generator_type,valence_generator_type,to;valenceflow=true)
        effOps = flow_Operators(binfo,HFobj,IMSRGobj,PandyaObj,Chan1b,Chan2bD,d9j,HOBs,dictMono,dict6j,Operators,to)
        if binfo.nuc.cZ != binfo.nuc.Z || binfo.nuc.cN != binfo.nuc.N
            println("1 Zerobody ",IMSRGobj.H.zerobody)
            getNormalOrderedO(binfo,HFobj,IMSRGobj.H,Chan1b,Chan2bD,dict6j,to;undo=true,OpeqH=true)
            println("2 undo =>  ",IMSRGobj.H.zerobody)
            set_sps_to_core!(binfo,HFobj)
            getNormalOrderedO(binfo,HFobj,IMSRGobj.H,Chan1b,Chan2bD,dict6j,to;OpeqH=true)            
            println("3 ms =>  ",IMSRGobj.H.zerobody)
            for Op in effOps
                set_sps_to_modelspace!(binfo,HFobj)
                getNormalOrderedO(binfo,HFobj,Op,Chan1b,Chan2bD,dict6j,to;undo=true,OpeqH=false)
                set_sps_to_core!(binfo,HFobj)
                getNormalOrderedO(binfo,HFobj,Op,Chan1b,Chan2bD,dict6j,to;OpeqH=false)
            end
        end
        write_vs_snt(binfo,HFobj,IMSRGobj,Operators,effOps,Chan1b,Chan2bD,valencespace)
    else
        flow_Operators(binfo,HFobj,IMSRGobj,PandyaObj,Chan1b,Chan2bD,d9j,HOBs,dictMono,dict6j,Operators,to)
    end
    return nothing
end

"""
    read_imsrg_parameter!(fn,IMSRGobj)
Function to overwrite IMSRGobj from the parameter file `fn`.
"""
function read_imsrg_parameter!(fn,IMSRGobj)
    include(fn)
    if @isdefined(smax); IMSRGobj.smax = smax; end
    if @isdefined(dsmax); IMSRGobj.dsmax = dsmax; end
    if @isdefined(maxnormOmega); IMSRGobj.maxnormOmega = maxnormOmega; end
    if @isdefined(eta_criterion); IMSRGobj.eta_criterion = eta_criterion; end
    if @isdefined(denominatorDelta); IMSRGobj.denominatorDelta = denominatorDelta;end
    println("parameters in $fn will be used.")
    return nothing
end

"""
    init_dictMonopole!(dictMonopole,Chan2b)

initialize `dictMonopole`
"""
function init_dictMonopole!(dictMonopole,Chan2b)
    for ch = 1:length(Chan2b)
        tbc = Chan2b[ch]
        kets = tbc.kets; J = tbc.J; Tz = tbc.Tz
        pnrank = 2 + div(Tz,2)
        tdict = dictMonopole[pnrank]
        for (idx,ket) in enumerate(kets)
            exists = get(tdict,ket,false)
            if exists == false
                tdict[ket] = valDictMonopole([0.0,0.0],[ [ch,idx,J] ])
            else
                push!(tdict[ket].vals,[ch,idx,J])
            end
        end
    end
    return nothing
end

"""
    set_dictMonopole!(dictMonopole,HFobj,H) 

To update dictMonopole pp/pn/nn under H(s=0)/IMSRG H(s)
"""
function set_dictMonopole!(dictMonopole,HFobj,H) 
    MS = HFobj.modelspace; sps = MS.sps
    for pnrank =1:3
        for key in keys(dictMonopole[pnrank])
            tmp = dictMonopole[pnrank][key]
            sqT0 = 1+delta(key[1],key[2])
            sqfac = ifelse(pnrank!=2,sqT0,1.0)
            monovec = tmp.monopole
            ja = sps[key[1]].j
            jb = sps[key[2]].j
            monovec .*= 0.0
            chidxs = tmp.vals
            for ch_idx_J in chidxs
                ch,idx,J = ch_idx_J
                NJ = 2*J + 1
                V = H[ch][idx,idx]
                monovec[1] += NJ * V * sqfac / ( (ja+1)*(jb+1)) 
            end 
        end
    end 
    return nothing
end

"""
    getNorm(O,p_sps,n_sps,Chan2b)

returns sqrt(norm1b^2 + norm2b^2)
"""
function getNorm(O,p_sps,n_sps,Chan2b)
    n1b = getNorm1b(O.onebody,p_sps,n_sps)
    n2b = getNorm2b(O.twobody,Chan2b)
    return sqrt(n1b^2+n2b^2)
end

"""
    getNorm1b(Mat1b,p_sps,n_sps,verbose=false)

returns 1bnorm of the given Operator
"""
function getNorm1b(Mat1b,p_sps,n_sps,verbose=false)
    tnorm =0.0
    for pn=1:2
        M = Mat1b[pn]
        dim = size(M)[1]
        sps = ifelse(pn==1,p_sps,n_sps)
        for i = 1:dim
            ji = sps[i].j
            for j=1:dim
                degen = ji + 1
                tmp = (M[i,j] * degen)^2
                tnorm += tmp
                if verbose && abs(M[i,j]) > 1.e-6
                    println("1bnorm single pn $pn i $i j $j tmp $tmp M[i,j] ",M[i,j],
                    "  degen  $degen")
                end
            end
        end
    end 
    return sqrt(tnorm)
end 

"""
    getNorm2b(Mat2b,Chan2b,verbose=false)

returns 2bnorm of the given Operator
"""
function getNorm2b(Mat2b,Chan2b,verbose=false)
    tnorm = 0.0
    for ch =1:length(Mat2b)
        J = Chan2b[ch].J
        tmp = norm(Mat2b[ch],2)
        tnorm += (tmp * (2*J+1))^2
        if verbose; println("ch $ch  tnorm $tnorm");end
    end 
    return sqrt(tnorm)
end 
"""
    calc_Eta_atan!(HFobj,IMSRGobj,Chan2b,dictMono,norms)

calc. ``\\eta(s)`` with atan generator
"""
function calc_Eta_atan!(HFobj,IMSRGobj,Chan2b,dictMono,norms)
    MS = HFobj.modelspace
    sps = MS.sps; holes = MS.holes; particles = MS.particles
    p_sps = MS.p_sps; n_sps = MS.n_sps
    key = zeros(Int64,2)
    ## one-body piece a->h i->v/p
    Hs = IMSRGobj.H; f = Hs.onebody;  Gamma = Hs.twobody
    Eta1b = IMSRGobj.eta.onebody
    Eta2b = IMSRGobj.eta.twobody
    Delta = IMSRGobj.denominatorDelta    
    for i=1:2
        eta1b = Eta1b[i]; eta1b .*= 0.0
    end
    for (a,oa) in enumerate(sps) # a->core
        if !oa.c;continue;end 
        idx_a = div(a,2) + a%2
        for (i,oi) in enumerate(sps) # i-> v/q (= not c)
            if a > i;continue;end
            if oa.tz != oi.tz;continue;end
            if oi.c; continue;end 
            pn = ifelse(oa.tz==-1,1,2)
            eta1b = Eta1b[pn]
            tf = f[pn]
            pnrank = ifelse(pn==1,1,3)
            dMono = dictMono[pnrank]
    
            idx_i = div(i,2) + i%2
            nume = 2 * tf[idx_a,idx_i]
            key[1] = a; key[2] = i
            if a > i; key[1] = i; key[2]= a; end
            mono_ai = dMono[key].monopole[1] 
            deno = tf[idx_i,idx_i] - tf[idx_a,idx_a] + (oi.occ-oa.occ)*mono_ai + Delta
            tmp = 0.5 * atan( nume / deno)
            eta1b[idx_i,idx_a] = tmp
            eta1b[idx_a,idx_i] = -tmp
        end
    end 
    ## two-body piece hh vs pp
    spaces = HFobj.modelspace.spaces
    cc = spaces.cc
    for ch in keys(cc)
        Gam = Gamma[ch]; tbc = Chan2b[ch]; Tz = tbc.Tz; kets = tbc.kets
        eta2b = Eta2b[ch]; eta2b .*= 0.0
        pnrank = div(Tz,2) + 2
        idxs_cc = cc[ch]
        #println("ch $ch nket ",length(kets))
        for ik in idxs_cc
            i,j = kets[ik]
            ni = sps[i].occ; nj = sps[j].occ
            for ib =1:length(kets)
                a,b = kets[ib]
                na = sps[a].occ; nb = sps[b].occ
                if sps[a].c || sps[b].c ; continue;end                
                nume = 2 * Gam[ib,ik]
                deno = Get2bDenominator(ch,pnrank,a,b,i,j,na,nb,ni,nj,f,Delta,dictMono,key)                             
                tmp = 0.5 * atan(nume / deno)                
                eta2b[ib,ik] = tmp
                if ib != ik; eta2b[ik,ib] = -tmp;end
            end
        end 
    end 
    p_sps = HFobj.modelspace.p_sps
    n_sps = HFobj.modelspace.n_sps
    norms[3] = getNorm1b(IMSRGobj.eta.onebody,p_sps,n_sps)
    norms[4] = getNorm2b(IMSRGobj.eta.twobody,Chan2b)    
    return nothing
end

"""
    check_order_Mkey(key,pnrank)

reorder key to be ```key[1] > key[2]```
"""
function check_order_Mkey(key,pnrank)
    i,j = key
    if i > j;
        key[1] = j; key[2]=i
    end
    return nothing
end 

"""
    Get2bDenominator(ch,pnrank,a,b,i,j,na,nb,ni,nj,f,Delta,dictMono,key;verbose=false)

`` f_{aa} +f_{bb} −f_{ii} −f_{jj} +G_{abij} +\\Delta``

with ``G_{abij} = \\Gamma_{abab} + \\Gamma_{ijij} - (\\Gamma_{aiai} + \\Gamma_{bjbj} + [a \\leftrightarrow b])``
"""
function Get2bDenominator(ch,pnrank,a,b,i,j,na,nb,ni,nj,f,Delta,dictMono,key;verbose=false)
    key[1] = a; key[2] = b; check_order_Mkey(key,pnrank); Vm_abab = dictMono[pnrank][key].monopole[1] # pp'-pp'
    key[1] = i; key[2] = j; check_order_Mkey(key,pnrank); Vm_ijij = dictMono[pnrank][key].monopole[1] # hh'-hh'
    pnrank_ai = 3 - (a%2 + i%2); pnrank_bj= 3 - (b%2 + j%2)
    key[1] = a; key[2] = i; check_order_Mkey(key,pnrank_ai); Vm_aiai = dictMono[pnrank_ai][key].monopole[1] # ph-ph
    key[1] = b; key[2] = j; check_order_Mkey(key,pnrank_bj); Vm_bjbj = dictMono[pnrank_bj][key].monopole[1] # p'h'-p'h'
    pnrank_aj = 3 - (a%2 + j%2); pnrank_bi = 3 - (b%2 + i%2)
    key[1] = i; key[2] = b; check_order_Mkey(key,pnrank_bi); Vm_bibi = dictMono[pnrank_bi][key].monopole[1] # p'h-p'h
    key[1] = a; key[2] = j; check_order_Mkey(key,pnrank_aj); Vm_ajaj = dictMono[pnrank_aj][key].monopole[1] # ph'-ph'
    idx_a = div(a,2) + a%2; idx_b = div(b,2) + b%2
    idx_i = div(i,2) + i%2; idx_j = div(j,2) + j%2
    pn_a = 2 - a%2; pn_b = 2 - b%2; pn_i = 2 - i%2; pn_j = 2 - j%2
    denoG  = f[pn_a][idx_a,idx_a] + f[pn_b][idx_b,idx_b] - f[pn_i][idx_i,idx_i] -f[pn_j][idx_j,idx_j]
    denoG += ( 1-na-nb ) * Vm_abab - (1-ni-nj) * Vm_ijij + Delta
    denoG += (na-ni) * Vm_aiai +  (nb-nj) * Vm_bjbj +  (nb-ni) * Vm_bibi +  (na-nj) * Vm_ajaj
    return denoG
end

"""
    print_flowstatus(istep,s,ncomm,norms,IMSRGobj)

print flowstatus s,E0,1b&2b norm for Omega, 1b&2b norm for Eta, Ncomm, nwritten
"""
function print_flowstatus(istep,s,ncomm,norms,IMSRGobj)
    if istep == 0
        println(" step:        s             E0   ||Omega_1||   ||Omega_2||     ||Eta_1||     ||Eta_2||    Ncomm.  nwritten")
    end        
    E0 = IMSRGobj.H.zerobody[1]
    println(@sprintf("%6i",istep),"  ",@sprintf("%7.3f",s), 
        @sprintf("%15.8f",E0),
        @sprintf("%14.6e",norms[1]),@sprintf("%14.6e",norms[2]),
        @sprintf("%14.6e",norms[3]),@sprintf("%14.6e",norms[4]),
        @sprintf("%10i",ncomm[1]),@sprintf("%10i",IMSRGobj.n_written_omega[1]))
    return nothing
end

"""
    make_PandyaKets(emax,HFobj)

To prepare "kets" for Pandya transformation.
For ordinary two-body channels, kets like ```|i,j=i;J=odd>``` where ```i={n,l,j,tz}``` are hindered, but necessary for Pandya transformation.
"""
function make_PandyaKets(emax,HFobj)
    Chan2b_Pandya = chan2b[ ]
    MS = HFobj.modelspace;sps = MS.sps
    dim1b = div(length(sps),2)
    Jmax = 2*emax+1
    nchan = nchanP = 0
    ns = Vector{Int64}[ ]
    ns_addinv = Vector{Int64}[ ]
    for Tz = -2:2:2 ## pp/pn/nn
        for prty = 1:-2:-1        
            for J = 0:Jmax
                kets = Vector{Int64}[ ]
                nhh = nph = 0
                for a = 1:2*dim1b
                    oa = sps[a]; la = oa.l; ja = oa.j; tza = oa.tz; na = oa.occ
                    for b=a:2*dim1b
                        ob = sps[b]; lb = ob.l; jb = ob.j; tzb = ob.tz; nb = ob.occ
                        aTz = abs(tza+tzb)
                        if aTz != abs(Tz); continue;end
                        if tri_check(ja//2,jb//2,J)==false;continue;end
                        tprty = (-1)^(la+lb)
                        if prty != tprty; continue;end
                        #if oa.occ + ob.occ == 2; nhh +=1;end                                             
                        #if oa.occ + ob.occ == 1; nph +=1;end
                        if (na != 0.0) && (nb !=0.0); nhh +=1;end 
                        if (na*nb ==0.0) && (na+nb!=0.0);nph +=1;end
                        push!(kets,[a,b])
                    end
                end
                nchan += 1
                nket = length(kets)
                if nket > 0 
                    if Tz == -2 # for addinv 
                        push!(ns_addinv,[nchan,nket,nhh,nph])
                    else
                        nchanP += 1
                        nket = length(kets)
                        push!(Chan2b_Pandya, chan2b(Tz,prty,J,kets))  
                        push!(ns,[nchan,nket,nhh,nph])
                        push!(ns_addinv,[nchan,nket,nhh,nph])
                    end
                end
            end
        end
    end
    ## sort Chan2bPandya by nKets_cc
    nkets = [ ns[i][2] for i=1:length(ns)]
    idxs = sortperm(nkets,rev=true)
    sortedChan2b = chan2b[]
    sorted_ns = Vector{Int64}[]
    dict_ch2ich = Dict{Int64,Int64}()
    for (ich,idx) in enumerate(idxs)
        push!(sortedChan2b,Chan2b_Pandya[idx])      
        push!(sorted_ns,ns[idx])
        ch = ns[idx][1]
        dict_ch2ich[ch] = ich
    end  
    return sorted_ns,ns_addinv,sortedChan2b,dict_ch2ich
end 

"""
constructor of utils for Pandya transformation and others
numbers_Pandya:[ch,nKet_cc,nhh,nph] for ich (channel index of Chan2b_Pandya) 
"""
function prep_PandyaLookup(binfo,HFobj,Chan1b,Chan2bD;rank_J=0,rank_T=0,parity=0,ofst=1000)   
    Chan2b = Chan2bD.Chan2b
    Jmax = binfo.emax * 2 + 1
    numbers_Pandya,numbers_forAddInv,Chan2b_Pandya,dict_ch2ich = make_PandyaKets(binfo.emax,HFobj)
    MS = HFobj.modelspace; sps = MS.sps

    ##prep occupation matrix for na*nb and nabar*nbbar
    sps = HFobj.modelspace.sps
    Mats_hh = Matrix{Float64}[ ]
    Mats_pp = Matrix{Float64}[ ]
    Mats_ph = Matrix{Float64}[ ]
    for ch = 1:length(Chan2b)
        tbc = Chan2b[ch];kets = tbc.kets
        ppidx = get(HFobj.modelspace.spaces.pp,ch,Int64[]); nhh = length(ppidx) 
        hhidx = get(HFobj.modelspace.spaces.hh,ch,Int64[]); nhh = length(hhidx) 
        phidx = get(HFobj.modelspace.spaces.ph,ch,Int64[]); nph = length(phidx) 
        Mat_pp = zeros(Float64,nhh,nhh);Mat_hh = zeros(Float64,nhh,nhh)
        Mat_ph =  zeros(Float64,nph,nph)
        for (i,idx) in enumerate(hhidx) #this is right! pp from hh!
            a,b = kets[idx]; na = sps[a].occ; nb = sps[b].occ
            Mat_pp[i,i] = (1-na)*(1-nb)
        end
        for (i,idx) in enumerate(hhidx)
            a,b = kets[idx]; na = sps[a].occ; nb = sps[b].occ
            Mat_hh[i,i] = na*nb
        end
        for (i,idx) in enumerate(phidx)
            a,b = kets[idx]; na = sps[a].occ; nb = sps[b].occ
            Mat_ph[i,i] = (1-na)*(1-nb) #+ na *nb
        end
        push!(Mats_hh,Mat_hh)
        push!(Mats_pp,Mat_pp)
        push!(Mats_ph,Mat_ph)
    end

    ## prep. XYbars for workers
    nchPandya = length(Chan2b_Pandya)
    Zbars = Matrix{Float64}[ ]
    phkets = Vector{Int64}[ ]
    PhaseMats = Matrix{Float64}[]   
    dimmax = 0
    dict_ich_idx_from_ketcc = [ Dict{Int64,Int64}( ) for ich=1:nchPandya]        
    for ich = 1:nchPandya
        ch,nKet_cc,nhh,npp = numbers_Pandya[ich]
        nph_kets = nhh + npp  
        push!(Zbars,zeros(Float64,nKet_cc,2*nKet_cc))
        dimmax = maximum([dimmax,nKet_cc,2*nph_kets])
        push!(phkets,Gethhph(Chan2b_Pandya[ich].kets,sps))
        push!(PhaseMats, zeros(Float64,nKet_cc,nKet_cc))
        tbc = Chan2b_Pandya[ich]; tkets = tbc.kets
        target = dict_ich_idx_from_ketcc[ich]
        for (idx,ket) in enumerate(tkets)
            tkey = ket[1] * ofst + ket[2]
            target[tkey] = idx
        end 
    end

    nthre = nthreads()
    XYbars =  [ [zeros(Float64,dimmax,dimmax),zeros(Float64,dimmax,dimmax)] for i=1:nthre]
    tMat = [ zeros(Float64,2*dimmax,2*dimmax) for i=1:nthre]
    keys6j = [zeros(Int64,5) for i=1:nthreads()]

    ### to prep. 122 util such as intermediate matrix
    util122 = prep122(HFobj,Chan1b,Chan2bD)
    return PandyaObject(numbers_Pandya,numbers_forAddInv,Chan2b_Pandya,phkets,                
                        dict_ich_idx_from_ketcc,XYbars,Zbars,PhaseMats,tMat,dict_ch2ich,keys6j,util122,Mats_hh,Mats_pp,Mats_ph)
end

function lookup_twobody(a,d,c,b,ja,jd,jc,jb,Jtar,lookup,O2b,Chan2b,key6j;verbose=false)
    phase = 1.0 
    key1 = @view key6j[1:2]; key1[1] = a; key1[2] = d
    key2 = @view key6j[3:4]; key2[1] = c; key2[2] = b
    if a > d; key1[1] = d; key1[2] = a; phase *= (-1)^(div(ja+jd,2)+Jtar+1); end
    if c > b; key2[1] = b; key2[2] = c; phase *= (-1)^(div(jb+jc,2)+Jtar+1); end
    ad_candidates = lookup[key1]
    cb_candidates = lookup[key2]
    #println("adcb $a $d $c $b  key1 $key1 key2 $key2 J $Jtar $ad_candidates")
    tbme = 0.0
    for ad_cand in ad_candidates
        ch,idx_ad = ad_cand        
        if Chan2b[ch].J != Jtar;continue;end
        for cb_cand in cb_candidates
            ch_cb,idx_cb = cb_cand
            if ch != ch_cb; continue;end
            ta,td = Chan2b[ch].kets[idx_ad]
            tc,tb = Chan2b[ch].kets[idx_cb]
            sqfac = ifelse(ta==td,sqrt(2.0),1.0) * ifelse(tc==tb,sqrt(2.0),1.0)
            tbme += O2b[ch][idx_ad,idx_cb] * phase * sqfac
            #println("hit@ch=$ch idx_ad $idx_ad ta td $ta $td  idx_cb $idx_cb tc tb $tc $tb")
        end
    end
    return tbme
end

function DoPandyaTransformation(O,Obar_ph,tbc_cc,Chan2bD,HFobj,PandyaObj,numbers_ch,dict6j,key6j,orientation="N")
    Chan2b = Chan2bD.Chan2b
    MS = HFobj.modelspace; sps = MS.sps
    J_cc = tbc_cc.J
    O2b = O.twobody
    herm = ifelse(O.hermite,1.0,-1.0)
    ch_cc,nKets_cc,nhh,nph = numbers_ch
    nph_kets = nhh + nph
    hit_ph = 0    
    tdict = dict6j[J_cc+1]
    Dict_chidx_from_ketJ = Chan2bD.dict_ch_idx_from_ket
    for ibra = 1:nKets_cc #ph_kets
        bra_cc = tbc_cc.kets[ibra]
        a,b = bra_cc
        if (sps[a].occ +sps[b].occ == 0.0);continue;end
        hit_ph += 1
        # to consider a > b case 
        for ab_case = 0:1 # 0-> |a,b> 1->|b,a>
            a = bra_cc[1 + ab_case]
            b = bra_cc[2 - ab_case]
            na = sps[a].occ; nb = sps[b].occ
            ja = sps[a].j; jb = sps[b].j
            na_nb_factor = na-nb
            bra_shift = 0
            if a == b && ab_case==1
                bra_shift=nph_kets
            end
            if (na != 0.0 && nb == 0.0)
                bra_shift = 0
            elseif (na ==0.0 && nb != 0.0)
                bra_shift = nph_kets
            else
                bra_shift = ifelse(ab_case==0,0,nph_kets)
            end
            for iket_cc = 1:nKets_cc
                ket_cc = tbc_cc.kets[iket_cc]
                c,d = ket_cc
                jc = sps[c].j; jd = sps[d].j
                Tz = sps[a].tz + sps[d].tz
                if (sps[a].tz + sps[d].tz) != (sps[b].tz + sps[c].tz);continue;end
                jmin = max(div(abs(ja-jd),2), div(abs(jc-jb),2))
                jmax = min(div(ja+jd,2),div(jc+jb,2))
                Obar = 0.0
                for J_std = jmin:jmax
                    if !tri_check(ja,jd,2*J_std);continue;end
                    if !tri_check(jb,jc,2*J_std);continue;end
                    if abs(Tz)==2 && J_std %2 ==1 && (b==c || a==d);continue;end
                    nkey = get_nkey_from_key6j(ja,jb,jc,jd,J_std)
                    sixj = tdict[nkey]
                    if abs(sixj) <= 1.e-8;continue;end
                    lookup = Dict_chidx_from_ketJ[1+div(Tz+2,2)][J_std+1]                    
                    tbme = lookup_twobody(a,d,c,b,ja,jd,jc,jb,J_std,lookup,O2b,Chan2b,key6j) 
                    Obar -= (2*J_std + 1) * sixj * tbme
                end
                if orientation == "N" # for Y
                    Obar_ph[hit_ph+bra_shift,iket_cc] = Obar
                elseif orientation =="T" # for X                    
                    Obar_ph[iket_cc,hit_ph+bra_shift] = herm * Obar * na_nb_factor
                else
                    println("err orientation=$orientation is not defined");exit()
                end
            end
        end
    end
    return nothing
end

function get_nkets_permutation(tbc)
    kets = tbc.kets
    Tz = tbc.Tz
    count = 0
    for ket in kets
        count +=1 
        if Tz != 0 #&& ket[1] != ket[2]
            count +=1
        end
    end
    return count
end


function need_check(tbc_bra_cc, tbc_ket_cc, li,ji,tzi, lj,jj,tzj, lk,jk,tzk, ll,jl,tzl)
    tf = false
    twoJ_bra_cc = tbc_bra_cc.J * 2
    twoJ_ket_cc = tbc_ket_cc.J * 2
    par_bra = tbc_bra_cc.prty
    par_ket = tbc_ket_cc.prty
    Tz_bra  = tbc_bra_cc.Tz
    Tz_ket  = tbc_ket_cc.Tz
    ## i-l & j-k
    j3min = abs(ji-jl); j3max = ji+jl
    j4min = abs(jk-jj); j4max = jk+jj
    if ( ((-1)^(li+ll) == par_bra) && 
         ((-1)^(lj+lk) == par_ket) && 
         (abs(tzi+tzl) == Tz_bra) && 
         (abs(tzj+tzk) == Tz_ket) && 
         j3min<=twoJ_bra_cc && twoJ_bra_cc<=j3max &&
         j4min<=twoJ_ket_cc && twoJ_ket_cc<=j4max )
        #if ch_bra_cc == ch_ket_cc == 1; println("hit @1");end
        return true
    end
    if (((-1)^(li+ll) == par_ket) && 
        ((-1)^(lj+lk) == par_bra) && 
        (abs(tzi+tzl) == Tz_ket)   && 
        (abs(tzj+tzk) == Tz_bra)   && 
        j3min<=twoJ_ket_cc && twoJ_ket_cc<=j3max &&
        j4min<=twoJ_bra_cc && twoJ_bra_cc<=j4max )
        #if ch_bra_cc == ch_ket_cc == 1; println("hit @2");end
       return true
   end
   ##i-k j-l
   j3min = abs(jj-jl); j3max = jj+jl
   j4min = abs(ji-jk); j4max = ji+jk
   if (((-1)^(lj+ll) == par_bra) && 
       ((-1)^(li+lk) == par_ket) && 
       (abs(tzj+tzl) == Tz_bra)   && 
       (abs(tzi+tzk) == Tz_ket)   && 
        j3min<=twoJ_bra_cc && twoJ_bra_cc<=j3max &&
        j4min<=twoJ_ket_cc && twoJ_ket_cc<=j4max )
        #if ch_bra_cc == ch_ket_cc == 1; println("hit @3");end
        return true
   end
   if (((-1)^(lj+ll) == par_ket) && 
       ((-1)^(li+lk) == par_bra) && 
       (abs(tzj+tzl) == Tz_ket)   && 
       (abs(tzi+tzk) == Tz_bra)   && 
        j3min<=twoJ_ket_cc && twoJ_ket_cc<=j3max &&
        j4min<=twoJ_bra_cc && twoJ_bra_cc<=j4max )
        #if ch_bra_cc == ch_ket_cc == 1; println("hit @4");end
        return true
   end
   return tf
end


"""
    adhoc_rewrite6jdict(emax,dict6j,ofst_unit=1000)

adhoc function to replace dict. for 6j-symbols:{ji,jj,J,jk,jl,J'}
key = [ji,jj,jk,jl,J'] -> new_key::Int64
new_dict is Vector for total J
"""
function adhoc_rewrite6jdict(emax,dict6j,ofst_unit=1000)
    Jmax = 2*emax+1 
    new_dict6j = [ Dict{Int64,Float64}() for tJ=0:Jmax]
    for totJ = 0:Jmax
        tdict  = dict6j[totJ+1]
        target = new_dict6j[totJ+1]        
        for tkey in keys(tdict)
            tval = tdict[tkey]
            nkey = get_nkey_from_key6j(tkey,ofst_unit)
            target[nkey] = tval
        end
    end
    return new_dict6j
end 
function  get_nkey_from_key6j(tkey,ofst_unit=1000)
    nkey = tkey[1] + tkey[2] * ofst_unit + tkey[3] * ofst_unit^2 + tkey[4] * ofst_unit^3 +  tkey[5] * ofst_unit^4
    return nkey
end
function get_nkey_from_key6j(ji,jj,jk,jl,Jp,ofst=1000)
    return ji + jj*ofst + jk * ofst^2+ jl* ofst^3 + Jp * ofst^4
end

function IMSRGflow(binfo,HFobj,IMSRGobj,PandyaObj,Chan1b,Chan2bD,dictMono,dict6j,
                   core_generator,valence_generator,to;valenceflow=false,debugmode=0,maxstep=2000,magnusmethod="split")        
    Chan2b = Chan2bD.Chan2b
    ncomm = IMSRGobj.Ncomm    
    s,ds = IMSRGobj.s
    smax = IMSRGobj.smax 
    maxnormOmega=IMSRGobj.maxnormOmega
    eta_criterion = IMSRGobj.eta_criterion
    norms = zeros(Float64,4)
    H0 = IMSRGobj.H0
    Hs = IMSRGobj.H
    Omega = IMSRGobj.Omega
    eta = IMSRGobj.eta
    func_Eta = ""
    if !valenceflow 
        if core_generator == "atan" 
            func_Eta = calc_Eta_atan!
        else
            println("core_generator=$core_generator is not supported now");exit()
        end
    else    
        println("\n\nStarting VS-IMSRG flow")
        aOp1_p_bOp2!(Hs,H0,1.0,0.0)
        aOp!(Omega,0.0)
        if valence_generator == "shell-model-atan"
            func_Eta = calc_Eta_smatan!
        else
            println("valence_generator=$valence_generator is not supported now");exit()
        end
    end 
    set_dictMonopole!(dictMono,HFobj,Hs.twobody)
    func_Eta(HFobj,IMSRGobj,Chan2b,dictMono,norms)
    tmpOp = deepcopy(IMSRGobj.Omega)
    gatherer = ifelse(magnusmethod=="huntergather",deepcopy(IMSRGobj.Omega),nothing)
    nOmega = deepcopy(IMSRGobj.Omega)
    Nested = deepcopy(IMSRGobj.Omega)
    istep = 0; print_flowstatus(istep,s,ncomm,norms,IMSRGobj) 
    @timeit to "IMSRG flow" for istep = 1:maxstep
        if sqrt(norms[3]^2+norms[4]^2) < eta_criterion;break;end
    	ds = min(ds,smax-s) #最後を考慮
    	s += ds
        IMSRGobj.s[1] = s; IMSRGobj.s[2] = ds
        #debug_print(debugmode,Hs,Omega,eta,"s")
    
        ## OmegaCheck
        if magnusmethod == "huntergather"
            GatherOmega(Omega,nOmega,gatherer,tmpOp,Nested,H0,Hs,
                        ncomm,norms,Chan1b,Chan2bD,HFobj,IMSRGobj,dictMono,dict6j,PandyaObj,maxnormOmega,to)
        elseif magnusmethod == "split" # not implemented
            NewOmega(binfo,Omega,nOmega,HFobj,IMSRGobj,Chan2bD)
        elseif magnusmethod == "" || magnusmethod == "NS"
            aOp1_p_bOp2!(nOmega,Omega,1.0,0.0)
        else
            println("magnusmethod=$magnusmethod is not supported!");exit()
        end
        # Eta * ds for Euler step (assuming Magnus expansion) 
        aOp!(eta,ds) 
        ## Calc Omega(s+ds) 
        BCH_Product(eta,Omega,nOmega,tmpOp,Nested,ncomm,norms,Chan1b,Chan2bD,HFobj,dictMono,dict6j,PandyaObj,to)

        ## IMSRGobj H(s+ds) updated 
        BCH_Transform(nOmega,H0,Hs,tmpOp,Nested,ncomm,norms,Chan1b,Chan2bD,HFobj,dictMono,dict6j,PandyaObj,to) 

        #debug_print(debugmode,Hs,nOmega,eta,"s+ds")
        set_dictMonopole!(dictMono,HFobj,IMSRGobj.H.twobody)
        func_Eta(HFobj,IMSRGobj,Chan2b,dictMono,norms)

        print_flowstatus(istep,s,ncomm,norms,IMSRGobj)
        if s >= smax;break;end
    end
    aOp1_p_bOp2!(nOmega,Omega,1.0,0.0)
    write_omega_bin(binfo,IMSRGobj.n_written_omega[1],nOmega)
    IMSRGobj.n_written_omega[1] += 1
    #show(to, allocations = true,compact = false);println("") 
    return nothing
end   
function NewOmega(binfo,Omega,nOmega,HFobj,IMSRGobj,Chan2bD)
    maxnormOmega = IMSRGobj.maxnormOmega
    MS = HFobj.modelspace; p_sps =MS.p_sps; n_sps=MS.n_sps
    Chan2b = Chan2bD.Chan2b
    dim1b = size(Omega.onebody[1])[1]
    if getNorm(nOmega,p_sps,n_sps,Chan2b) >= maxnormOmega  
        write_omega_bin(binfo,IMSRGobj.n_written_omega[1],nOmega)
        IMSRGobj.n_written_omega[1] += 1
        aOp!(Omega,0.0)
        H0 = IMSRGobj.H0
        Hs = IMSRGobj.H
        aOp1_p_bOp2!(Hs,H0,1.0,0.0)
    else
        aOp1_p_bOp2!(nOmega,Omega,1.0,0.0)
    end
    return nothing
end

"""
    write_omega_bin(binfo,n_written,Omega)
to write binary file of Operator matrix elements, when spliting the flow
"""
function write_omega_bin(binfo,n_written,Omega)
    if isdir("flowOmega") 
        if n_written==0
            rm.(glob("flowOmega/*bin"))
        end
    else
        run(`mkdir flowOmega`)
    end    
    pid = getpid()
    nw = n_written + 1
    fname = "flowOmega/Omega_$pid"*binfo.nuc.cnuc*"_$nw.bin"
    dim1b = size(Omega.onebody[1])[1]
    nch = length(Omega.twobody)
    dims = Int64[]
    for ch = 1:nch
        dim = size(Omega.twobody[ch])[1]
        push!(dims,dim)
    end
    io = open(fname,"w")
    write(io,Int64(dim1b))
    for i =1:dim1b
        for j=1:dim1b
           write(io,Omega.onebody[1][i,j])
        end
    end
    for i =1:dim1b
        for j=1:dim1b
            write(io,Omega.onebody[2][i,j])
        end
    end
    write(io,Int64(nch))
    write(io,dims)
    for ch = 1:nch
        dim = dims[ch]
        O2b = Omega.twobody[ch]
        for i =1:dim
            for j=1:dim
                write(io,O2b[i,j])
            end
        end
    end 
    close(io)
    return nothing
end

"""
    read_omega_bin!(nw,Op,verbose=false)

read written Omega file and update ```Op::Operator```
"""
function read_omega_bin!(binfo,nw,Op,verbose=false)
    pid = getpid()
    fn = "flowOmega/Omega_$pid"*binfo.nuc.cnuc*"_$nw.bin"
    io = open(fn,"r")
    dim1b = read(io,Int64)
    for i = 1:dim1b
        for j=1:dim1b
            Op.onebody[1][i,j] = read(io,Float64)
        end
    end
    for i = 1:dim1b
        for j=1:dim1b
            Op.onebody[2][i,j] = read(io,Float64)
        end
    end
    nch = read(io,Int64)
    dims = [read(io,Int64) for ch=1:nch]
    for ch = 1:nch
        dim2b = dims[ch]
        O2b = Op.twobody[ch]
        for i=1:dim2b
            for j=1:dim2b
                O2b[i,j] = read(io,Float64)
            end
        end
    end 
    close(io)
    return nothing
end

"""
    GatherOmega(Omega,nOmega,gatherer,tmpOp,Nested,H0,Hs,ncomm,norms,Chan1b,Chan2bD,HFobj,IMSRGobj,dictMono,dict6j,PandyaObj,maxnormOmega,to)

This may not be used now.
"""
function GatherOmega(Omega,nOmega,gatherer,tmpOp,Nested,H0,Hs,ncomm,norms,Chan1b,Chan2bD,HFobj,IMSRGobj,dictMono,dict6j,PandyaObj,maxnormOmega,to)
    MS = HFobj.modelspace; p_sps =MS.p_sps; n_sps=MS.n_sps
    Chan2b = Chan2bD.Chan2b
    maxnormOmega = IMSRGobj.maxnormOmega
    if getNorm(Omega,p_sps,n_sps,Chan2b) <= maxnormOmega
        aOp1_p_bOp2!(nOmega,Omega,1.0,0.0)
        return false
    end    
    aOp!(gatherer,0.0)
    BCH_Product(nOmega,Omega,gatherer,tmpOp,Nested,ncomm,norms,Chan1b,Chan2bD,HFobj,dictMono,dict6j,PandyaObj,to)
    BCH_Transform(gatherer,H0,Hs,tmpOp,Nested,ncomm,norms,Chan1b,Chan2bD,HFobj,dictMono,dict6j,PandyaObj,to)        
    set_dictMonopole!(dictMono,HFobj,IMSRGobj.H.twobody)
    aOp1_p_bOp2!(gatherer,Omega,1.0,0.0)
    return true    
end

function debug_print(debugmode,Hs,Omega,eta,label;ch=1)
    if debugmode == 2
        label *= ";ch1"
        print_vec("Omega($label ,ch=$ch) \t",Omega.twobody[ch][1,:])
        print_vec("Eta($label ,ch=$ch) \t",eta.twobody[ch][1,:])
    elseif debugmode == 1
        println("proton 1b($label):") 
        for i=1:size(Hs.onebody[1])[1]
            print_vec("",Hs.onebody[1][i,:])
        end
    end
    return nothing
end

"""
    Gethhph(kets,sps)

get idxs for hh/ph kets
"""
function Gethhph(kets,sps)
    idx_hhph = Int64[ ]
    for (idx,ket) in enumerate(kets)
        p,q = ket 
        occs = sps[p].occ + sps[q].occ 
        #if occs == 2 || occs == 1; push!(idx_hhph,idx);end
        if occs != 0.0; push!(idx_hhph,idx);end # for fractional
    end
    return idx_hhph
end

"""
    init_IMSRGobject(HFobj;smax=500.0,dsmax=0.5,maxnormOmega=0.25,tol=1.e-6,eta_criterion=1.e-6,denominatorDelta=0.0)

Constructor for IMSRGobject
- `H0::Operator` for starting point of BCH product
- `H::Operator` Hamiltonian ``H(s)``
- `s::Vector{Float}` current ``s`` and ``ds``
- `smax::Float` maximum ``s``
- `dsmax::Float` maximum ``ds``
- `maxnormOmega::Float` maximum ||Omega||
- `eta::Operator` generator of IMSRG flow (antihermite Operator)
- `Omega::Operator` generator of IMSRG flow (antihermite Operator) 
- `eta_criterion::Float` ||eta|| to check convergence
- `denominatorDelta::Float64` parameter for multi-major shell decoupling
- `n_written_omega::Int` # of written Omega by splitting to solve IMSRGflow
- `Ncomm::Vector{Int}` # of commutator evaluated during IMSRG flow
"""
function init_IMSRGobject(HFobj;smax=500.0,dsmax=0.5,maxnormOmega=0.25,eta_criterion=1.e-6,denominatorDelta=0.0,filename="optional_parameters.jl")
    tf = isfile(filename)
    ds = 0.5
    s = [0.0,ds] 
    E0 = HFobj.E0
    fp = copy(HFobj.H.onebody[1])
    fn = copy(HFobj.H.onebody[2])
    Gamma = deepcopy(HFobj.H.twobody)
    Hs = Operator([E0],[fp,fn],Gamma,true,false)
    H0 = deepcopy(Hs)
    eta = Operator([0.0],[copy(0.0 .*fp),copy(0.0 .*fn)],0.0 .*deepcopy(Gamma),false,true)
    Omega = Operator([0.0],[copy(0.0 .*fp),copy(0.0 .*fn)], 0.0 .*deepcopy(Gamma),false,true)
    n_written_omega = [0]
    Ncomm =[0] 

    IMSRGobj = IMSRGobject(H0,Hs,s,smax,dsmax,maxnormOmega,eta,Omega,eta_criterion,denominatorDelta,n_written_omega,Ncomm)
    if !tf
        println("Since $filename is not found, the default parameters will be used.")
        return IMSRGobj
    else
        read_imsrg_parameter!(filename,IMSRGobj)
        return IMSRGobj
    end
end

"""
    flow_Operators(binfo,HFobj,IMSRGobj,PandyaObj,Chan1b,Chan2bD,dict_9j,HOBs,dictMono,dict6j,Operators,to)

consistent IMSRG flow of scaler operators (Rp2) using written Omega
"""
function flow_Operators(binfo,HFobj,IMSRGobj,PandyaObj,Chan1b,Chan2bD,dict_9j,HOBs,dictMono,dict6j,Operators,to)
    effOperators = Operator[ ]
    for c_op in Operators
        if c_op=="Rp2"
            println("Operator:$c_op")
            eOp = eval_rch_imsrg(binfo,Chan1b,Chan2bD,HFobj,IMSRGobj,PandyaObj,dict_9j,HOBs,dictMono,dict6j,to)
            push!(effOperators,eOp)
        else
            println("IMSRG flow of Operator=$c_op is not supported now ")
            continue
        end
    end
    return effOperators
end

"""
    update_core_in_sps!(binfo,HFobj) 

Function to specify hole/core for sps. This will will be used for target normal ordering
"""
function update_core_in_sps!(binfo,HFobj)    
    cZ = binfo.nuc.cZ
    cN = binfo.nuc.cN
    sps = HFobj.modelspace.sps
    p_sps = HFobj.modelspace.p_sps
    n_sps = HFobj.modelspace.n_sps
    ### to specify "core" occupation
    dim1b = length(p_sps)
    occ_p = zeros(Float64,dim1b,dim1b)
    occ_n = zeros(Float64,dim1b,dim1b)
    pconfs_core = naive_filling(p_sps,cZ,binfo.emax)
    nconfs_core = naive_filling(n_sps,cN,binfo.emax)    
    ini_occ!(pconfs_core,occ_p,nconfs_core,occ_n)
    for i = 1:dim1b
        pTF = p_sps[i].occ==1.0 
        nTF = n_sps[i].occ==1.0
        p_sps[i].c = sps[2*(i-1)+1].c = ifelse(pTF,true,false)
        n_sps[i].c = sps[2*i].c = ifelse(nTF,true,false)
    end
    return nothing
end

"""
    set_sps_to_core!(binfo,HFobj)

modify ```p_sps.occ, n_sps.occ``` by the "core" nucleus
"""
function set_sps_to_core!(binfo,HFobj)
    cZ = binfo.nuc.cZ
    cN = binfo.nuc.cN
    sps = HFobj.modelspace.sps
    p_sps = HFobj.modelspace.p_sps
    n_sps = HFobj.modelspace.n_sps
    dim1b = length(p_sps)
    occ_p = zeros(Float64,dim1b,dim1b)
    occ_n = zeros(Float64,dim1b,dim1b)
    pconfs_core = naive_filling(p_sps,cZ,binfo.emax)
    nconfs_core = naive_filling(n_sps,cN,binfo.emax)    
    ini_occ!(pconfs_core,occ_p,nconfs_core,occ_n)
    for i = 1:dim1b
        cp = ifelse(occ_p[i,i]==1.0,1,0)
        cn = ifelse(occ_n[i,i]==1.0,1,0)
        p_sps[i].occ = sps[2*(i-1)+1].occ = cp
        n_sps[i].occ = sps[2*i].occ = cn
    end
    return nothing
end

"""
    set_sps_to_modelspace!(binfo,HFobj)

modify occupation by specified model space
"""
function set_sps_to_modelspace!(binfo,HFobj)
    cZ = binfo.nuc.Z
    cN = binfo.nuc.N
    sps = HFobj.modelspace.sps
    p_sps = HFobj.modelspace.p_sps
    n_sps = HFobj.modelspace.n_sps
    dim1b = length(p_sps)
    occ_p = zeros(Float64,dim1b,dim1b)
    occ_n = zeros(Float64,dim1b,dim1b)
    pconfs_core = naive_filling(p_sps,cZ,binfo.emax)
    nconfs_core = naive_filling(n_sps,cN,binfo.emax)    
    ini_occ!(pconfs_core,occ_p,nconfs_core,occ_n)
    for i = 1:dim1b
        cp = ifelse(occ_p[i,i]==1.0,1,0)
        cn = ifelse(occ_n[i,i]==1.0,1,0)
        p_sps[i].occ = sps[2*(i-1)+1].occ = cp
        n_sps[i].occ = sps[2*i].occ = cn
    end
    return nothing
end
