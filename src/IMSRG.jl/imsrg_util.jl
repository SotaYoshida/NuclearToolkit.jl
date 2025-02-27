"""
    imsrg_main(binfo::basedat,Chan1b,Chan2bD,HFobj,dictMono,dWS,valencespace,Operators,to; core_generator_type="atan",valence_generator_type="shell-model-atan",denominatorDelta=0.0)
    
# Arguments
- `binfo::basedat` struct basedat(nuc::nuclei,sntf::String,hw::Float,emax::Int)
- `Chan1b::chan1b` struct for one-body stuffs
- `Chan2bD::chan2bD` struct for two-body stuffs (e.g., dict to get idx from JPT)
- `HFobj::HamiltonianNormalOrdered` struct HNO, which includes info. of HF solution (HF energy, occupation, f,Gamma,...)
- `dictMono::Dict` dictionary to get Vmonopole
- `dWS` dictionary of preallocated wigner-symbols
- `valencespace` to specify valence space  
- `Operators::Vector{String}` non-Hamiltonian operators
- `to` TimerOutput object to measure runtime&memory allocations

# Optional Arguments
- `delete_Ops` if true, delete Operators with current pid after IMSRGflow
- `core_generator_type` only the "atan" is available
- `valence_generator_type` only the "shell-model-atan" is available
- `denominatorDelta::Float` denominator Delta, which is needed for multi-major shell decoupling
- `debugmode=0`: 0: no debug, 1: debug, 2: debug with more info
- `restart_from_files`: files to be read for restart 1st one is for IMSRG and 2nd one is for VSIMSRG
"""
function imsrg_main(binfo::basedat,Chan1b::chan1b,Chan2bD::chan2bD,HFobj::HamiltonianNormalOrdered,dictsnt,dWS,valencespace,Operators,MatOp,to;
                    delete_Ops=false,core_generator_type="atan",valence_generator_type="shell-model-atan",fn_params="optional_parameters.jl",debugmode=0,Hsample=0,restart_from_files=String[])
    vsIMSRG = ifelse(valencespace!=[]&&valencespace!="",true,false)
    if binfo.nuc.corenuc == "" && vsIMSRG; println("core (hole) for VS-IMSRG is not specified."); return nothing;end 
    dictMono = deepcopy(dictsnt.dictMonopole)
    update_core_in_sps!(binfo,HFobj)
    valencesps = check_valence_space(HFobj,valencespace)
    update_vsspace_chs!(HFobj,valencesps,Chan2bD.Chan2b)

    Chan2b = Chan2bD.Chan2b
    init_dictMonopole!(dictMono,Chan2b)
    IMSRGobj = init_IMSRGobject(HFobj,fn_params)
    PandyaObj = prep_PandyaLookup(binfo,HFobj,Chan1b,Chan2bD)
    if length(restart_from_files) >= 1; IMSRGobj.smax = IMSRGobj.dsmax;end
    d6j_defbyrun = Dict{UInt64,Float64}()
    IMSRGflow(binfo,HFobj,IMSRGobj,PandyaObj,Chan1b,Chan2bD,dictMono,d6j_defbyrun,core_generator_type,valence_generator_type,to;debugmode=debugmode,Hsample=Hsample,restart_from_files=restart_from_files)
    IMSRGobj.ExpectationValues["E0"] = IMSRGobj.H.zerobody[1]

    if vsIMSRG && length(restart_from_files) == 0
        IMSRGflow(binfo,HFobj,IMSRGobj,PandyaObj,Chan1b,Chan2bD,dictMono,d6j_defbyrun,core_generator_type,valence_generator_type,to;valenceflow=true,Hsample=Hsample)
        effOps = flow_Operators(binfo,HFobj,IMSRGobj,PandyaObj,Chan1b,Chan2bD,dWS,d6j_defbyrun,dictMono,Operators,MatOp,restart_from_files,to)
        if binfo.nuc.cZ != binfo.nuc.Z || binfo.nuc.cN != binfo.nuc.N
            getNormalOrderedO(HFobj,IMSRGobj.H,Chan1b,Chan2bD,to;undo=true,OpeqH=true)
            set_sps_to_core!(binfo,HFobj)
            getNormalOrderedO(HFobj,IMSRGobj.H,Chan1b,Chan2bD,to;OpeqH=true)            
            for Op in effOps
                set_sps_to_modelspace!(binfo,HFobj)
                getNormalOrderedO(HFobj,Op,Chan1b,Chan2bD,to;undo=true,OpeqH=false)
                set_sps_to_core!(binfo,HFobj)
                getNormalOrderedO(HFobj,Op,Chan1b,Chan2bD,to;OpeqH=false)
            end
        end
        write_vs_snt(binfo,HFobj,IMSRGobj,Operators,Chan1b,Chan2bD,valencespace;effOps=effOps)
    else        
        flow_Operators(binfo,HFobj,IMSRGobj,PandyaObj,Chan1b,Chan2bD,dWS,d6j_defbyrun,dictMono,Operators,MatOp,restart_from_files,to)
    end
    
    if length(restart_from_files) >= 1
        if binfo.nuc.cZ != binfo.nuc.Z || binfo.nuc.cN != binfo.nuc.N
            println("normal ordering...")
            getNormalOrderedO(HFobj,IMSRGobj.H,Chan1b,Chan2bD,to;undo=true,OpeqH=true)
            println("onebody@1 ", IMSRGobj.H.onebody[1][1,:])
            set_sps_to_core!(binfo,HFobj)
            getNormalOrderedO(HFobj,IMSRGobj.H,Chan1b,Chan2bD,to;OpeqH=true)            
            println("onebody@2 ", IMSRGobj.H.onebody[1][1,:])
            for Op in effOps
                set_sps_to_modelspace!(binfo,HFobj)
                getNormalOrderedO(HFobj,Op,Chan1b,Chan2bD,to;undo=true,OpeqH=false)
                set_sps_to_core!(binfo,HFobj)
                getNormalOrderedO(HFobj,Op,Chan1b,Chan2bD,to;OpeqH=false)
            end
        end
        #write_vs_snt(binfo,HFobj,IMSRGobj,Operators,Chan1b,Chan2bD,valencespace)
    end
    if delete_Ops
        pid = getpid()
        for f in glob("flowOmega/Omega_$(pid)*.bin")
            rm(f)
        end
    end
    return IMSRGobj
end

"""
    read_imsrg_parameter!(fn,IMSRGobj)
Function to overwrite IMSRGobj from the parameter file `fn`.
"""
function read_imsrg_parameter!(fn::String,IMSRGobj::IMSRGobject)
    include(fn)
    if @isdefined(smax); IMSRGobj.smax = smax; end
    if @isdefined(dsmax); IMSRGobj.dsmax = dsmax; end
    if @isdefined(maxnormOmega); IMSRGobj.maxnormOmega = maxnormOmega; end
    if @isdefined(eta_criterion); IMSRGobj.eta_criterion = eta_criterion; end
    if @isdefined(denominatorDelta); IMSRGobj.denominatorDelta = denominatorDelta;end
    if @isdefined(magnusmethod); IMSRGobj.magnusmethod = magnusmethod;end
    IMSRGobj.s[2] = IMSRGobj.dsmax
    println("parameters in $fn will be used.")
    return nothing
end

"""
    init_dictMonopole!(dictMonopole,Chan2b)

initialize `dictMonopole`
"""
function init_dictMonopole!(dictMonopole,Chan2b)
    for ch in eachindex(Chan2b)
        tbc = Chan2b[ch]
        kets = tbc.kets; J = tbc.J; Tz = tbc.Tz
        pnrank = 2 + div(Tz,2)
        tdict = dictMonopole[pnrank]
        for (idx,ket) in enumerate(kets)
            tkey = zeros(Int64,2)
            tkey .= ket
            if !haskey(tdict,tkey)
                tdict[tkey] = valDictMonopole([0.0,0.0],[ [ch,idx,J] ])
            else   
                push!(tdict[tkey].vals,[ch,idx,J])
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
            sqT0 = 1.0+delta(key[1],key[2])
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
function getNorm(O::Operator,p_sps,n_sps,Chan2b)
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
function getNorm2b(Mat2b::Vector{Matrix{Float64}},Chan2b::Vector{chan2b},verbose=false)
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
    calc_Eta_atan!(HFobj::HamiltonianNormalOrdered,IMSRGobj::IMSRGobject,Chan2b::Vector{chan2b},dictMono,norms)

calc. ``\\eta(s)`` with atan generator
"""
function calc_Eta_atan!(HFobj::HamiltonianNormalOrdered,IMSRGobj::IMSRGobject,Chan2b::Vector{chan2b},dictMono,norms)
    MS = HFobj.modelspace
    sps = MS.sps; p_sps = MS.p_sps; n_sps = MS.n_sps
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
        if !oa.c[1];continue;end 
        idx_a = div(a,2) + a%2
        for (i,oi) in enumerate(sps) # i-> v/q (= not c)
            if a > i;continue;end
            if oa.tz != oi.tz;continue;end
            if oi.c[1]; continue;end 
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
            deno = tf[idx_i,idx_i] - tf[idx_a,idx_a] + (oi.occ[1]-oa.occ[1])*mono_ai + Delta
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
            ni = sps[i].occ[1]; nj = sps[j].occ[1]
            for ib in eachindex(kets)
                a,b = kets[ib]
                na = sps[a].occ[1]; nb = sps[b].occ[1]
                if sps[a].c[1] || sps[b].c[1] ; continue;end                
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

`` f_{aa}+f_{bb}-f_{ii}-f_{jj}+G_{abij} +\\Delta``

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

Function to print flowstatus s,E0,1b&2b norm for Omega, 1b&2b norm for Eta, Ncomm, nwritten
"""
function print_flowstatus(istep,s,ncomm,norms,IMSRGobj,Chan2b,verbose=false)
    if istep == 0
        if verbose
            println(" # of 2b-ch ",length(IMSRGobj.Omega.twobody))
            dim_mat = dim_vec = 0
            for ch = 1:length(Chan2b)
                tbc = Chan2b[ch]
                d = tbc.nkets
                dim_mat += d^2
                dim_vec += div(d*(d+1),2)
            end
            println("dim mat $dim_mat fvec $dim_vec")
        end
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
For ordinary two-body channels, kets like `|i,j=i;J=odd>` with ``={n,l,j,tz}` are hindered, but necessary for Pandya transformation.
"""
function make_PandyaKets(emax::Int,HFobj::HamiltonianNormalOrdered,Chan2b)
    Chan2b_Pandya = chan2b[ ]
    MS = HFobj.modelspace;sps = MS.sps
    dim1b = div(length(sps),2)
    Jmax = 2*emax+1
    nchan = nchanP = 0
    ns = NTuple{4,Int64}[ ]
    ns_addinv = NTuple{4,Int64}[ ]
    
    for Tz = -2:2:2 ## pp/pn/nn
        for prty = 1:-2:-1        
            for J = 0:Jmax
                kets = NTuple{2,Int64}[ ]                
                nhh = nph = 0
                for a = 1:2*dim1b
                    oa = sps[a]; la = oa.l; ja = oa.j; tza = oa.tz; na = oa.occ[1]
                    for b=a:2*dim1b
                        ob = sps[b]; lb = ob.l; jb = ob.j; tzb = ob.tz; nb = ob.occ[1]
                        aTz = abs(tza+tzb)
                        if aTz != abs(Tz); continue;end
                        if tri_check(ja//2,jb//2,J)==false;continue;end
                        tprty = (-1)^(la+lb)
                        if prty != tprty; continue;end
                        if (na != 0.0) && (nb !=0.0); nhh +=1;end 
                        if (na*nb ==0.0) && (na+nb!=0.0);nph +=1;end
                        push!(kets,(a,b))
                    end
                end
                nchan += 1
                nket = length(kets)                
                if nket > 0
                    if Tz == -2 # for addinv
                        push!(ns_addinv,(nchan,nket,nhh,nph))                   
                    else
                        nchanP += 1                        
                        push!(Chan2b_Pandya, chan2b(Tz,prty,J,kets,nket))  
                        push!(ns,(nchan,nket,nhh,nph))
                        push!(ns_addinv,(nchan,nket,nhh,nph))
                    end
                end
            end
        end
    end
    ## sort Chan2bPandya by nKets_cc
    nkets = [ ns[i][2] for i in eachindex(ns)]
    idxs = sortperm(nkets,rev=true)
    sortedChan2b = chan2b[]
    sorted_ns = NTuple{4,Int64}[]
    dict_ch2ich = zeros(Int64,length(Chan2b))
    for (ich,idx) in enumerate(idxs)
        push!(sortedChan2b,Chan2b_Pandya[idx])      
        push!(sorted_ns,ns[idx])
        ch = ns[idx][1]
        dict_ch2ich[ch] = ich
    end      
    return sorted_ns,ns_addinv,sortedChan2b,dict_ch2ich
end 

function prep_nab_bar_matrices(HFobj::HamiltonianNormalOrdered,Chan2bD::chan2bD)
    Chan2b = Chan2bD.Chan2b
    MS = HFobj.modelspace; sps = MS.sps
    nch = length(Chan2b)
    mats_nab = Vector{Float64}[ ]
    mats_nab_bar = Vector{Float64}[ ]
    for ch = 1:nch
        tkets = Chan2b[ch].kets
        nket = length(tkets)
        mat_nab = zeros(Float64,nket)
        mat_nab_bar = zeros(Float64,nket)
        @inbounds for idx_ab = 1:nket
            a,b = tkets[idx_ab]
            na = sps[a].occ[1]
            nb = sps[b].occ[1]
            nab = 1.0 * (na * nb)
            nab_bar = 1.0 * (1-na)*(1-nb)
            mat_nab[idx_ab] = nab
            mat_nab_bar[idx_ab] = nab_bar
        end
        push!(mats_nab,mat_nab)
        push!(mats_nab_bar,mat_nab_bar)
    end 
    return mats_nab,mats_nab_bar
    # mats_* are now vectors though...
end

"""
    prep_PandyaLookup(binfo::basedat,HFobj::HamiltonianNormalOrdered,Chan1b::chan1b,Chan2bD::chan2bD;rank_J=0,rank_T=0,parity=0,ofst=1000)   

constructor of utils for Pandya transformation and others
numbers_Pandya:[ch,nKet_cc,nhh,nph] for ich (channel index of Chan2b_Pandya) 
"""
function prep_PandyaLookup(binfo::basedat,HFobj::HamiltonianNormalOrdered,Chan1b::chan1b,Chan2bD::chan2bD;rank_J=0,rank_T=0,parity=0,ofst=1000)   
    Chan2b = Chan2bD.Chan2b    
    numbers_Pandya,numbers_forAddInv,Chan2b_Pandya,dict_ch2ich = make_PandyaKets(binfo.emax,HFobj,Chan2b)
    MS = HFobj.modelspace; sps = MS.sps

    ##prep occupation matrix for na*nb and nabar*nbbar
    sps = HFobj.modelspace.sps
    Mats_hh = Matrix{Float64}[ ]
    Mats_pp = Matrix{Float64}[ ]
    Mats_ph = Matrix{Float64}[ ]
    for ch in eachindex(Chan2b)
        tbc = Chan2b[ch];kets = tbc.kets
        ppidx = get(HFobj.modelspace.spaces.pp,ch,Int64[]); nhh = length(ppidx) 
        hhidx = get(HFobj.modelspace.spaces.hh,ch,Int64[]); nhh = length(hhidx) 
        phidx = get(HFobj.modelspace.spaces.ph,ch,Int64[]); nph = length(phidx) 
        Mat_pp = zeros(Float64,nhh,nhh);Mat_hh = zeros(Float64,nhh,nhh)
        Mat_ph =  zeros(Float64,nph,nph)
        for (i,idx) in enumerate(hhidx) #this is right! pp from hh!
            a,b = kets[idx]; na = sps[a].occ[1]; nb = sps[b].occ[1]
            Mat_pp[i,i] = (1-na)*(1-nb)
        end
        for (i,idx) in enumerate(hhidx)
            a,b = kets[idx]; na = sps[a].occ[1]; nb = sps[b].occ[1]
            Mat_hh[i,i] = na*nb
        end
        for (i,idx) in enumerate(phidx)
            a,b = kets[idx]; na = sps[a].occ[1]; nb = sps[b].occ[1]
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
    dict_ich_idx_from_ketcc = [ Dict{UInt64,Int64}( ) for ich=1:nchPandya]
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
            tkey = get_nkey2_u(ket[1],ket[2])
            target[tkey] = idx
        end 
    end
    nthre = Threads.maxthreadid()
    dim1b = length(MS.p_sps)
    copy_1bmat = [ zeros(Float64,dim1b,dim1b) for i=1:2*nthre]
    XYbars =  [ [zeros(Float64,dimmax,dimmax),zeros(Float64,dimmax,dimmax)] for i=1:nthre]
    tMat = [ zeros(Float64,2*dimmax,2*dimmax) for i=1:nthre]
    keys6j = [zeros(Int64,5) for i=1:Threads.maxthreadid()]
    ## for 221
    mats_nab,mats_nab_bar = prep_nab_bar_matrices(HFobj,Chan2bD)
    ### to prep. 122 util such as intermediate matrix
    util122 = prep122(HFobj,Chan1b,Chan2bD)



    Pandya3Channels = NTuple{3,Int64}[]
    idxdict_nth = Dict{UInt64,Int64}()
    hermite = true
    for ich in eachindex(numbers_forAddInv)
        ch = numbers_forAddInv[ich][1]
        tbc = Chan2b[ch]; nKets = tbc.nkets
        tkets = tbc.kets
        for ibra in eachindex(tkets)
            ketmin = ifelse(hermite,ibra,ibra+1)
            for iket = ketmin:nKets
                push!(Pandya3Channels, (ch,ibra,iket))
                idxdict_nth[get_nkey3_u(ch,ibra,iket)] = length(Pandya3Channels)
            end
        end
    end   
    vec_flat = zeros(Float64, length(Pandya3Channels))
    return PandyaObject(numbers_Pandya,numbers_forAddInv,Chan2b_Pandya,phkets,copy_1bmat,mats_nab,mats_nab_bar,
                        dict_ich_idx_from_ketcc,XYbars,Zbars,PhaseMats,tMat,dict_ch2ich,keys6j,util122,Mats_hh,Mats_pp,Mats_ph,
                        vec_flat,idxdict_nth,Pandya3Channels)
end

function DoPandyaTransformation(O::Operator,Obar_ph,tbc_cc,Chan2bD,HFobj,nph_kets,d6j_lj,to,orientation="N";def_mode=false)
    MS = HFobj.modelspace; sps = MS.sps
    J_cc = tbc_cc.J
    tkets = tbc_cc.kets
    O2b = O.twobody
    herm = ifelse(O.hermite[1],1.0,-1.0)
    ch_dict = Chan2bD.dict_ch_idx_from_ket
    hit_ph = 0 
    @inbounds for ibra in eachindex(tkets)
        bra_cc = tkets[ibra]
        a,b = bra_cc
        if (sps[a].occ[1] +sps[b].occ[1] == 0.0);continue;end
        hit_ph += 1
        # to consider a > b case 
        for ab_case = 0:1 # 0-> |a,b> 1->|b,a>
            a = bra_cc[1 + ab_case]
            b = bra_cc[2 - ab_case]
            na = sps[a].occ[1]; ja = sps[a].j
            nb = sps[b].occ[1]; jb = sps[b].j
            na_nb_factor = na-nb
            bra_shift = 0
            if a == b && ab_case==1
                bra_shift = nph_kets
            end
            if (na != 0.0 && nb == 0.0)
                bra_shift = 0
            elseif (na ==0.0 && nb != 0.0)
                bra_shift = nph_kets
            else
                bra_shift = ifelse(ab_case==0,0,nph_kets)
            end
            for iket_cc in eachindex(tkets)
                c,d = tkets[iket_cc]
                Tz = sps[a].tz + sps[d].tz
                if Tz != (sps[b].tz + sps[c].tz);continue;end               
                tdict = ch_dict[div(Tz+2,2)+1]
                jc = sps[c].j; jd = sps[d].j
                jmin = max(div(abs(ja-jd),2), div(abs(jc-jb),2))
                jmax = min(div(ja+jd,2),div(jc+jb,2))
                tf_skip = (Tz!=0 && (b==c || a==d))
                jmin = ifelse(tf_skip&&jmin%2==1,jmin+1,jmin)
                jmax = ifelse(tf_skip&&jmax%2==1,jmax-1,jmax)
                if jmin > jmax; continue;end                
                jstep = ifelse(tf_skip,2,1)
                Obar = inner_sum_Pandya(jmin,jstep,jmax,a,b,c,d,ja,jb,jc,jd,J_cc,Tz,tdict,O2b,d6j_lj;def_mode=def_mode)

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

function flip_phase(o_a,o_b,j_a,j_b,J)::Float64
    if o_a <= o_b
        return 1.0
    else
        jflip_odd = (div(j_a+j_b,2)+J+1)%2 == 1
        return ifelse(jflip_odd,-1.0,1.0)
    end
end

function inner_sum_Pandya(jmin,jstep,jmax,a,b,c,d,ja,jb,jc,jd,J_cc,Tz,subdict,
                          O2b::Vector{Matrix{Float64}},d6j_lj;def_mode=false)
    Obar = 0.0
    J_cc2 = J_cc*2
    sqfac = ifelse(a==d,sqrt(2.0),1.0) * ifelse(c==b,sqrt(2.0),1.0) # equiv to t-counter part
    flip_ad = ifelse(a>d,true,false); flip_cb = ifelse(c>b,true,false)
    for J_std = jmin:jstep:jmax
        sixj = call_d6j_defbyrun(ja,jb,J_cc2,jc,jd,J_std*2,d6j_lj;def_mode=def_mode)
        if abs(sixj) < 1.e-8;continue;end        
        tbme = lookup_twobody(a,d,c,b,flip_ad,flip_cb,J_std,Tz,subdict,O2b)         
        tbme *= sixj 
        tbme *= flip_phase(a,d,ja,jd,J_std) * flip_phase(c,b,jc,jb,J_std)
        Obar -= hahat(J_std) * tbme 
    end
    return Obar * sqfac
end

function lookup_twobody(a::Int64,d::Int64,c::Int64,b::Int64,flip_ad::Bool,flip_cb::Bool,J::Int64,Tz::Int64,
                        lookup::Dict{UInt64,NTuple{2,Int64}},O2b::Vector{Matrix{Float64}})::Float64
    nkey = UInt64(0)
    if flip_ad 
        nkey = get_nkey3_ketJ(d,a,J)
    else
        nkey = get_nkey3_ketJ(a,d,J)
    end
    ch,idx_ad = lookup[nkey]
    if flip_cb 
        nkey = get_nkey3_ketJ(b,c,J)
    else
        nkey = get_nkey3_ketJ(c,b,J)
    end
    idx_cb = lookup[nkey][2]
    return O2b[ch][idx_ad,idx_cb] 
end

function sort_d6j_lj(d6j)
    ret = Dict{UInt64,Float64}()

    allkeys = sort(collect(keys(d6j)))
    for tkey in allkeys
        ret[tkey] = d6j[tkey]
    end
    ## key check
    count = 0
    for (ith,tkey) in enumerate(keys(d6j))
        skey = allkeys[ith]
        println("count $count tkey $tkey skey $skey")
        count += 1
        if count == 10;break;end            
    end
    return ret
end

function IMSRGflow(binfo::basedat,HFobj::HamiltonianNormalOrdered,IMSRGobj::IMSRGobject,PandyaObj::PandyaObject,Chan1b::chan1b,Chan2bD,dictMono,d6j_lj,
                   core_generator,valence_generator,to;valenceflow=false,debugmode=0,maxstep=2000,Hsample=0,
                   restart_from_files=String[])     
    Chan2b = Chan2bD.Chan2b
    ncomm = IMSRGobj.Ncomm
    s,ds = IMSRGobj.s
    smax = IMSRGobj.smax 
    dsmax = IMSRGobj.dsmax 
    maxnormOmega=IMSRGobj.maxnormOmega
    magnusmethod = IMSRGobj.magnusmethod
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

    tmpOp = deepcopy(IMSRGobj.Omega)
    gatherer = ifelse(magnusmethod=="huntergather",deepcopy(IMSRGobj.Omega),nothing)
    nOmega = deepcopy(IMSRGobj.Omega)
    Nested = deepcopy(IMSRGobj.Omega)

    ## To define by run: d6j_lj
    if !valenceflow
        def_d6j_lj_by_run!(eta,Omega,HFobj,Chan2bD,d6j_lj,PandyaObj,to)
        #d6j_lj = sort_d6j_lj(d6j_lj)
        println("def-by-run d6j_lj done! ", length(keys(d6j_lj)))
    end
    
    dict_idx_op_to_flatvec, dict_idx_flatvec_to_op, dict_if_idx_for_hdf5 = get_non0omega_idxs(HFobj,nOmega,Chan2b;debugmode=debugmode)   
    # writing out util. dict. to see correspondance between idx and ket
    if debugmode > 0
        write_dicts_idx_ket(HFobj, binfo, Chan2b, dict_idx_op_to_flatvec, dict_idx_flatvec_to_op, dict_if_idx_for_hdf5)
    end

    fvec = get_fvec_from_Op(s, nOmega, dict_idx_op_to_flatvec, dict_idx_flatvec_to_op)
    set_dictMonopole!(dictMono,HFobj,Hs.twobody)
    func_Eta(HFobj,IMSRGobj,Chan2b,dictMono,norms)

    istep = 0; print_flowstatus(istep,s,ncomm,norms,IMSRGobj,Chan2b)
    @timeit to "IMSRG flow" for istep = 1:maxstep
        ds = min(ds,smax-s)
        s += ds
        IMSRGobj.s[1] = s; IMSRGobj.s[2] = ds

        ## OmegaCheck
        if magnusmethod == "huntergather" 
            GatherOmega(Omega,nOmega,gatherer,tmpOp,Nested,H0,Hs,
                        ncomm,norms,Chan1b,Chan2bD,HFobj,IMSRGobj,dictMono,d6j_lj,PandyaObj,maxnormOmega,to)
        elseif magnusmethod == "" || magnusmethod == "split"
            NewOmega(s, Omega, nOmega, binfo, Chan2bD, HFobj, IMSRGobj)
        elseif magnusmethod == "no-split" || magnusmethod == "NS"
            aOp1_p_bOp2!(nOmega,Omega,1.0,0.0)
        else
            println("magnusmethod=$magnusmethod is not supported!");exit()
        end
        # Eta * ds for Euler step (assuming Magnus expansion) 
        aOp!(eta,ds) 

        ## Calc. Omega(s+ds) 
        BCH_Product(eta,Omega,nOmega,tmpOp,Nested,ncomm,norms,Chan1b,Chan2bD,HFobj,dictMono,d6j_lj,PandyaObj,to)

        ## IMSRGobj.H updated to H(s+ds)
        BCH_Transform(nOmega,H0,Hs,tmpOp,Nested,ncomm,norms,Chan1b,Chan2bD,HFobj,dictMono,d6j_lj,PandyaObj,to) 

        set_dictMonopole!(dictMono,HFobj,IMSRGobj.H.twobody)
        func_Eta(HFobj,IMSRGobj,Chan2b,dictMono,norms)

        # remnant for IMSRG-Net sampling
        if Hsample > 0 
            cond_1 = Hsample == 1 && ( (s <= 20.00  || s == 30.0  || s == 50.0 ) || valenceflow) # for DMD
            cond_2 = Hsample == 2 && ( (15.0 <= s <= 20.00  || s == 30.0  || s == 50.0 ) || valenceflow) # for IMSRG-Net
            if cond_1 || cond_2
                Nested = deepcopy(IMSRGobj.Omega)
                gather_omega_sofar_write(Hsample,istep, s, fvec, Omega, nOmega, tmpOp, binfo, Chan1b, Chan2bD, HFobj, IMSRGobj, dictMono, d6j_lj, PandyaObj,to,dict_idx_op_to_flatvec, dict_idx_flatvec_to_op,dict_if_idx_for_hdf5)
            end
        end

        print_flowstatus(istep,s,ncomm,norms,IMSRGobj,Chan2b)
        if sqrt(norms[3]^2+norms[4]^2) < eta_criterion || s >= smax
            aOp1_p_bOp2!(nOmega,Omega,1.0,0.0)            
            write_omega_bin(binfo,Chan2b,IMSRGobj.n_written_omega[1],nOmega,s,IMSRGobj.H.zerobody[1])
            gather_omega_sofar_write(Hsample,istep+1, s, fvec, Omega, nOmega, tmpOp, binfo, Chan1b, Chan2bD, HFobj, IMSRGobj, dictMono, d6j_lj, PandyaObj,to,dict_idx_op_to_flatvec, dict_idx_flatvec_to_op,dict_if_idx_for_hdf5)
            #svd_Op_twobody(s,Omega,Chan2b;verbose=true,max_rank=20)
            IMSRGobj.n_written_omega[1] += 1
            break
        end
    end

    if length(restart_from_files) >= 1
        s = 0.0
        fvec = zeros(Float64,1)
        emutype = ""
        fn1 = restart_from_files[1][1]
        if occursin("annomega",fn1)
            emutype = ";IMSRGNet"
        elseif occursin("dmd",fn1)
            emutype = ";DMD"
        end
        for (nth,fn) in enumerate(restart_from_files[1])
            inttype = 
            aOp!(Omega,0.0)
            s = parse(Float64,split(split(split(fn,"_")[end],"s")[end],".h5")[1])
            if occursin("dmdvec",fn) # make fvec from dmdvec
                dmdvec = read_dmdvec_hdf5(fn)
                write_fvec_hdf5(binfo,dmdvec,dict_if_idx_for_hdf5,s,IMSRGobj.H.zerobody[1];label="omega_dmd")
                pid = getpid()
                fn = "flowOmega/omega_dmd_vec_$pid"*binfo.nuc.cnuc*"_s"*strip(@sprintf("%6.2f",s))*".h5"
            end
            fvec = read_fvec_hdf5(fn)
            update_Op_with_fvec!(fvec,Omega,dict_idx_flatvec_to_op)
            BCH_Transform(Omega,HFobj.H,tmpOp,nOmega,Nested,ncomm,norms,Chan1b,Chan2bD,HFobj,dictMono,d6j_lj,PandyaObj,to)
            println("IMSRG from file $fn")
            println("En(s=",strip(@sprintf("%8.2f",s)),"$emutype) = ",tmpOp.zerobody[1])    
            if length(restart_from_files) == 2
                if length(restart_from_files[2]) == 0; continue; end 
                fn = restart_from_files[2][nth]
                fvec = read_fvec_hdf5(fn)
                update_Op_with_fvec!(fvec,Omega,dict_idx_flatvec_to_op)
                println("Onebody from fvec ", Omega.onebody[1][1,:])
                BCH_Transform(Omega,Hcopy,tmpOp,nOmega,Nested,ncomm,norms,Chan1b,Chan2bD,HFobj,dictMono,d6j_lj,PandyaObj,to)
                println("p1b=> ", tmpOp.onebody[1][1,:])  
                println("VSIMSRG from file $fn \nEn(s=",strip(@sprintf("%8.2f",s)),") = ",tmpOp.zerobody[1])
                aOp1_p_bOp2!(tmpOp,Hs,1.0,0.0)
            end
        end
    end
    return nothing
end 

function normHs(Hs::Operator)
    norm2b = sqrt(Hs.zerobody[1]^2)
    for ch = 1:length(Hs.twobody)
        norm2b += norm(Hs.twobody[ch],2)^2
    end
    norm2b = sqrt(norm2b)
    return sqrt(norm(Hs.onebody[1],2)^2 + norm(Hs.onebody[2],2)^2 + norm2b^2)
end

function gather_omega_sofar_write(Hsample, istep, s, fvec, oOmega, nOmega, tmpOp, binfo, Chan1b, Chan2bD,
                                  HFobj, IMSRGobj, dictMono, d6j_lj, PandyaObj,to,
                                  dict_idx_op_to_flatvec, dict_idx_flatvec_to_op,dict_if_idx_for_hdf5;debug_mode=true)
    if Hsample==0; return nothing;end
    magnusmethod = IMSRGobj.magnusmethod
    splitting = ifelse((magnusmethod!="" && magnusmethod!="split"),false,true)
    fvec .*= 0.0
    ncom = zeros(Int64,3)
    norms = zeros(Float64,3)
    tildeO = deepcopy(nOmega); aOp!(tildeO,0.0)
    tmpOp = deepcopy(tildeO)
    tmpOp2 = deepcopy(tildeO)
    Nested = deepcopy(tildeO)
    tOmega = deepcopy(tildeO)
    if splitting 
        nw = IMSRGobj.n_written_omega[1]    
        for i = 1:nw
            aOp!(tOmega,0.0)
            read_omega_bin!(binfo,Chan2bD.Chan2b,i,tOmega)
            BCH_Product(tOmega,tildeO,tmpOp,tmpOp2,Nested,ncom,norms,Chan1b,Chan2bD,HFobj,dictMono,d6j_lj,PandyaObj,to)
            aOp1_p_bOp2!(tmpOp,tildeO,1.0,0.0)
        end
        BCH_Product(nOmega,tildeO,tmpOp,tmpOp2,Nested,ncom,norms,Chan1b,Chan2bD,HFobj,dictMono,d6j_lj,PandyaObj,to)
        aOp1_p_bOp2!(tmpOp,tildeO,1.0,0.0)
    else
        aOp1_p_bOp2!(nOmega,tildeO,1.0,0.0)
    end
    get_fvec_from_Op!(s, fvec, tildeO, dict_idx_op_to_flatvec, dict_idx_flatvec_to_op)
    write_fvec_hdf5(binfo,fvec,dict_if_idx_for_hdf5,s,IMSRGobj.H.zerobody[1];label="omega")  
    if Hsample == 2
        get_fvec_from_Op!(s, fvec, IMSRGobj.eta, dict_idx_op_to_flatvec, dict_idx_flatvec_to_op)
        write_fvec_hdf5(binfo,fvec,dict_if_idx_for_hdf5,s,IMSRGobj.H.zerobody[1];label="eta")
    end
    return nothing 
end

function NewOmega(s, Omega, nOmega, binfo, Chan2bD, HFobj, IMSRGobj)       
    maxnormOmega = IMSRGobj.maxnormOmega
    MS = HFobj.modelspace; p_sps =MS.p_sps; n_sps=MS.n_sps
    Chan2b = Chan2bD.Chan2b
    if getNorm(nOmega,p_sps,n_sps,Chan2b) >= maxnormOmega  
        write_omega_bin(binfo,Chan2b,IMSRGobj.n_written_omega[1],nOmega,s,IMSRGobj.H.zerobody[1])
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
Function to write temporary binary files of Operator matrix elements, when spliting the flow.
"""
function write_omega_bin(binfo::basedat,Chan2b::Vector{chan2b},n_written::Int,Omega::Operator,s::Float64,E0::Float64;Oplabel="Omega")
    if !isdir("flowOmega")  && n_written==0
        mkdir("flowOmega")
    end    
    pid = getpid()
    nw = n_written + 1
    fname = "flowOmega/$(Oplabel)_$pid"*binfo.nuc.cnuc*"_$nw.bin"
    dim1b = size(Omega.onebody[1])[1]
    nch = length(Omega.twobody)
    dims = Int64[]
    JPTzs = Vector{Int64}[ ]
    for ch = 1:nch
        dim = size(Omega.twobody[ch])[1]
        push!(dims,dim)
        J = Chan2b[ch].J
        P = Chan2b[ch].prty
        Tz = Chan2b[ch].Tz
        push!(JPTzs,[J,P,Tz])
    end
    io = open(fname,"w")
    write(io,s)
    write(io,E0)
    write(io,Int64(dim1b))
    count1b = 0
    for i =1:dim1b
        for j=1:dim1b
           p1b = Omega.onebody[1][i,j]
           write(io,p1b)
           if i <= j && p1b != 0.0
                #println("i $i j $j p1b $p1b")
                count1b += 1
           end
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
        for n = 1:3
            write(io,Int64(JPTzs[ch][n]))
        end
    end
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

    # dim_pp = dim_nn = dim_pn = 0
    # for ch = 1:nch
    #     dim = dims[ch]
    #     J = Chan2b[ch].J
    #     P = Chan2b[ch].prty
    #     Tz = Chan2b[ch].Tz
    #     O2b = Omega.twobody[ch]
    #     #println("ch $ch JPT $J $P $Tz")
    #     if Tz == -2; dim_pp += dim*(dim+1)÷2;end
    #     if Tz == 0; dim_pn += dim*(dim+1)÷2;end
    #     if Tz == 2; dim_nn += dim*(dim+1)÷2;end
    #     if Tz != -2; continue;end
    #     for chp = 1:nch
    #         Jp = Chan2b[chp].J
    #         Pp = Chan2b[chp].prty
    #         Tzp = Chan2b[chp].Tz
    #         if J != Jp; continue; end
    #         if P != Pp; continue; end
    #         if Tz != - Tzp;continue;end
    #         O2bp = Omega.twobody[chp]
    #         #println("s $s ch $ch J $J P $P O-O' ",norm(O2b-O2bp,2))
    #     end
    # end 
    #println("dim: $(2*count1b+dim_pp+dim_pn+dim_nn) 1b $count1b 2bsum $(dim_pp+dim_pn+dim_nn) pp $dim_pp pn $dim_pn nn $dim_nn")
    return nothing
end



"""
    read_omega_bin!(nw,Op,verbose=false)

read written Omega file and update ```Op::Operator```
"""
function read_omega_bin!(binfo::basedat,Chan2b::Vector{chan2b},nw::Int,Op::Operator,verbose=false;Oplabel="Omega",fn="")
    aOp!(Op,0.0)
    pid = getpid()
    if fn == ""
        fn = "flowOmega/$(Oplabel)_$pid"*binfo.nuc.cnuc*"_$nw.bin"
    end
    io = open(fn,"r")
    s = read(io,Float64)
    E = read(io,Float64)
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
    JPTz = [ [read(io,Int64),read(io,Int64),read(io,Int64)] for ch=1:nch]
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
    GatherOmega(Omega,nOmega,gatherer,tmpOp,Nested,H0,Hs,ncomm,norms,Chan1b,Chan2bD,HFobj,IMSRGobj,dictMono,d6j_lj,PandyaObj,maxnormOmega,to)

This may not be used now.
"""
function GatherOmega(Omega::Op,nOmega::Op,gatherer::Op,tmpOp::Op,Nested::Op,
                     H0,Hs,ncomm,norms,Chan1b::chan1b,Chan2bD,HFobj,IMSRGobj,dictMono,d6j_lj,PandyaObj,maxnormOmega,to) where Op<:Operator
    MS = HFobj.modelspace; p_sps =MS.p_sps; n_sps=MS.n_sps
    Chan2b = Chan2bD.Chan2b
    maxnormOmega = IMSRGobj.maxnormOmega
    if getNorm(Omega,p_sps,n_sps,Chan2b) <= maxnormOmega
        aOp1_p_bOp2!(nOmega,Omega,1.0,0.0)
        return false
    end    
    aOp!(gatherer,0.0)
    BCH_Product(nOmega,Omega,gatherer,tmpOp,Nested,ncomm,norms,Chan1b,Chan2bD,HFobj,dictMono,d6j_lj,PandyaObj,to)
    BCH_Transform(gatherer,H0,Hs,tmpOp,Nested,ncomm,norms,Chan1b,Chan2bD,HFobj,dictMono,d6j_lj,PandyaObj,to)        
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
        occs = sps[p].occ[1] + sps[q].occ[1]
        if occs != 0.0; push!(idx_hhph,idx);end
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
function init_IMSRGobject(HFobj,filename;smax=500.0,dsmax=0.5,maxnormOmega=0.25,eta_criterion=1.e-6,denominatorDelta=0.0)
    tf = isfile(filename)
    ds = min(dsmax,1.0)
    s = [0.0,ds] 
    E0 = HFobj.E0
    fp = copy(HFobj.H.onebody[1])
    fn = copy(HFobj.H.onebody[2])
    Gamma = deepcopy(HFobj.H.twobody)
    Hs = Operator([E0],[fp,fn],Gamma,[true],[false])
    H0 = deepcopy(Hs)
    eta = Operator([0.0],[copy(0.0 .*fp),copy(0.0 .*fn)],0.0 .*deepcopy(Gamma),[false],[true])
    Omega = Operator([0.0],[copy(0.0 .*fp),copy(0.0 .*fn)], 0.0 .*deepcopy(Gamma),[false],[true])
    n_written_omega = [0]
    Ncomm =[0] 
    magnusmethod = ""
    ExpectationValues = Dict{String,Float64}()
    IMSRGobj = IMSRGobject(H0,Hs,s,smax,dsmax,maxnormOmega,magnusmethod,eta,Omega,eta_criterion,denominatorDelta,n_written_omega,Ncomm,ExpectationValues)
    if !tf
        println("Since $filename is not found, the default parameters will be used.")
        return IMSRGobj
    else
        println("$filename is used for IMSRG")
        read_imsrg_parameter!(filename,IMSRGobj)
        return IMSRGobj
    end
end

"""
    flow_Operators(binfo,HFobj,IMSRGobj,PandyaObj,Chan1b,Chan2bD,d6j_lj,dictMono,Operators,to)

consistent IMSRG flow of scaler operators (Rp2) using written Omega
"""
function flow_Operators(binfo,HFobj,IMSRGobj,PandyaObj,Chan1b::chan1b,Chan2bD,dWS,d6j_defbyrun,dictMono,Operators,MatOp,restart_from_files,to)
    effOperators = Operator[ ]
    for c_op in Operators
        if c_op=="Rp2"
            println("Operator:$c_op")
            eOp = eval_rch_imsrg(binfo,Chan1b,Chan2bD,HFobj,IMSRGobj,PandyaObj,dWS,d6j_defbyrun,dictMono,MatOp,restart_from_files,to)
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
        pTF = p_sps[i].occ[1]==1.0 
        nTF = n_sps[i].occ[1]==1.0
        p_sps[i].c[1] = sps[2*(i-1)+1].c[1] = ifelse(pTF,true,false)
        n_sps[i].c[1] = sps[2*i].c[1] = ifelse(nTF,true,false)
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
        cp = ifelse(occ_p[i,i]==1.0,1.0,0.0)
        cn = ifelse(occ_n[i,i]==1.0,1.0,0.0)
        p_sps[i].occ[1] = sps[2*(i-1)+1].occ[1] = cp
        n_sps[i].occ[1] = sps[2*i].occ[1] = cn
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
        cp = ifelse(occ_p[i,i]==1.0,1.0,0.0)
        cn = ifelse(occ_n[i,i]==1.0,1.0,0.0)
        p_sps[i].occ[1] = sps[2*(i-1)+1].occ[1] = cp
        n_sps[i].occ[1] = sps[2*i].occ[1] = cn
    end
    return nothing
end
