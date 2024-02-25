"""
    hf_main(nucs,sntf,hw,emax;verbose=false,Operators=String[],is_show=false,doIMSRG=false,valencespace=[],corenuc="",ref="nucl")

Main API to carry out HF/HFMBPT or IMSRG calculation from snt file
# Arguments
- `nucs::Vector{String}` target nuclei
- `sntf` path to input interaction file
- `hw` hbar omega
- `emax_calc` emax for HF/IMSRG

# Optional Arguments
- `verbose=false` to see detailed stdout for HF
- `Operators=String[]` target observables other than Hamiltonian
- `is_show=false` to show TimerOutput log (summary of run time and memory allocation)
- `doIMSRG=false` to carry out IMSRG/VSIMSRG calculation 
- `valencespace=""` to spacify the valence space (e.g., "sd-shell" or ["sd-shell"], [[0,1,1,-1],[0,1,3,-1], [0,1,1,1],[0,1,3,1]]), if this is not empty, it tries to perform VS-IMSRG calculations
- `corenuc=""` core nucleus, example=> "He4"
- `ref="nucl"` to specify target reference state, "core" or "nucl" is supported.
- `return_obj=false` to return `hfdata` or `imsrgdata` object from this function
- `oupfn=""` to specify output file (writing stdout) name
- `fn_params="optional_parameters.jl"` to specify the name of file to read optional parameters
- `debugmode=0` to specify debug mode (0: no debug, 1: debug, 2: debug with more details)
- `Hsample=false` to specify whether to sample IMSRG omega and eta operators for emulators in hdf5 format 
- `restart_from_files=String[]` to specify the files to restart IMSRG flow from (e.g., ["ann_omega_vec_s20.h5"]). If this has two elements, the first one is for IMSRG and the other one is for VS-IMSRG.
"""
function hf_main(nucs,sntf,hw,emax_calc;verbose=false,Operators=String[],is_show=false,
                doIMSRG=false,delete_Ops=false,valencespace="",corenuc="",ref="nucl",return_obj=false,oupfn="",
                fn_params="optional_parameters.jl",debugmode=0,Hsample=false,restart_from_files=String[],
                e1max_file=0, e2max_file=0, e3max_file=0, e3max=0, fn_3nf="")
    to = TimerOutput()
    BLAS.set_num_threads(1)
    if debugmode != 0
        println("BLAS.get_config() ",BLAS.get_config())
        println("BLAS.get_num_threads() $(BLAS.get_num_threads()) nthreads $(Base.Threads.nthreads())")
    end
    @assert isfile(sntf) "sntf:$sntf is not found!"
    
    io = select_io(false,"",nucs;use_stdout=true,fn=oupfn)
    chiEFTparams = init_chiEFTparams(;io=nothing)
    HFdata = prepHFdata(nucs,ref,["E"],corenuc)
    @timeit to "prep dWS2n" begin
        no_need_9j_HOB = ifelse("Rp2" in Operators, false, true)
        dWS = prep_dWS2n(chiEFTparams,to;emax_calc=max(emax_calc,e1max_file),no_need_9j_HOB=no_need_9j_HOB)
    end
    @timeit to "read NNint" begin
        TF = occursin(".bin",sntf)
        tfunc = ifelse(TF,readsnt_bin,readsnt)     
        nuc = def_nuc(nucs[1],ref,corenuc)
        binfo = basedat(nuc,sntf,hw,emax_calc,ref)
        sps,dicts1b,dicts = tfunc(sntf,binfo,to)
        A=nuc.A
        BetaCM = chiEFTparams.BetaCM
        Hamil,dictsnt,Chan1b,Chan2bD,Gamma,maxnpq = store_1b2b(sps,dicts1b,dicts,binfo)
        HCM = InitOp(Chan1b,Chan2bD.Chan2b)
        TCM = InitOp(Chan1b,Chan2bD.Chan2b)
        VCM = InitOp(Chan1b,Chan2bD.Chan2b)
        E0cm = 1.5 * BetaCM * hw
        if BetaCM !=0.0
            Calculate_RCM(binfo,Chan1b,Chan2bD.Chan2b,sps,VCM,dWS,to;non0_ij=false)           
            fac_HCM = 0.5 * BetaCM * Mm * hw^2 / (hc^2)
            aOp!(VCM,fac_HCM)            
            aOp1_p_bOp2!(VCM,HCM,1.0,0.0)
            CalculateTCM!(TCM,binfo,Chan1b,Chan2bD.Chan2b,sps)
            update_dicts_withHCM!(HCM,Chan2bD,dicts)
        end 
        MatOp = Matrix{Float64}[]
        if "Rp2" in Operators
            MatOp = [ zeros(Float64,maxnpq,maxnpq) for i=1:2*nthreads()]
        end
    end

    # reading three-body interaction in me3j.gz format generated with NuHamil
    if fn_3nf != "" && !isfile(fn_3nf)
        println("fn_3nf:$fn_3nf is not found!")
        exit()
    end
    @timeit to "read 3NF" Object_3NF = main_read_me3j(fn_3nf, emax_calc, e1max_file, e2max_file, e3max, e3max_file, sps, dWS, to)

    Aold = A
    for (i,tnuc) in enumerate(nucs)
        nuc = def_nuc(tnuc,ref,corenuc); A=nuc.A   
        binfo = basedat(nuc,sntf,hw,emax_calc,ref)
        print(io,"target: $tnuc Ref. => Z=$(nuc.Z) N=$(nuc.N) ")
        if BetaCM !=0.0 && Aold != A
            difA_RCM(VCM,Aold,A)
            aOp1_p_bOp2!(VCM,HCM,1.0,0.0)
            difA_TCM(TCM,Aold,A)
            update_dicts_withHCM!(HCM,Chan2bD,dicts)
        end 
        recalc_v!(A,dicts)
        Hamil,dictsnt,Chan1b,Chan2bD,Gamma,maxnpq = store_1b2b(sps,dicts1b,dicts,binfo)    
        if i > 1
            update_1b!(binfo,sps,Hamil)
            update_2b!(binfo,sps,Hamil,dictsnt.dictTBMEs,Chan2bD,dicts)
        end
        addHCM1b!(Hamil,HCM,A)
        addHCM1b!(Hamil,TCM)
        @timeit to "HF" begin 
            HFobj = hf_iteration(binfo,HFdata[i],sps,Hamil,dictsnt.dictTBMEs,Chan1b,Chan2bD,Gamma,maxnpq,dWS,to;verbose=verbose,Object_3NF=Object_3NF,io=io,E0cm=E0cm) 
        end
        if doIMSRG
           IMSRGobj = imsrg_main(binfo,Chan1b,Chan2bD,HFobj,dictsnt,dWS,valencespace,Operators,MatOp,to;delete_Ops=delete_Ops,fn_params=fn_params,debugmode=debugmode,Hsample=Hsample,restart_from_files=restart_from_files)
           if return_obj; return IMSRGobj;end
        else
            if "Rp2" in Operators
                @timeit to "Rp2" begin 
                    Op_Rp2 = InitOp(Chan1b,Chan2bD.Chan2b)
                    eval_rch_hfmbpt(binfo,Chan1b,Chan2bD,HFobj,Op_Rp2,dWS,MatOp,to;io=io)
                end
            end
            if return_obj; return HFobj;end
        end
        Aold = A
    end
    show_TimerOutput_results(to;io=io,tf=is_show)
    return true
end

"""
    hf_main_mem(chiEFTobj,nucs,dict_TM,dWS,to;verbose=false,Operators=String[],valencespace="",corenuc="",ref="core")    
"without I/O" version of `hf_main`
"""
function hf_main_mem(chiEFTobj::ChiralEFTobject,nucs,dict_TM,dWS,HFdata,to;verbose=false,Operators=String[],valencespace="",corenuc="",ref="core",io=stdout) 
    emax = chiEFTobj.params.emax
    hw = chiEFTobj.params.hw
    sntf = chiEFTobj.params.fn_tbme   
    nuc = def_nuc(nucs[1],ref,corenuc)
    Z = nuc.Z; N=nuc.N; A=nuc.A 
    binfo = basedat(nuc,sntf,hw,emax,ref)
    sps,p_sps,n_sps = def_sps(emax)
    new_sps,dicts1b = make_sps_and_dict_isnt2ims(p_sps,n_sps,emax)
    dicts = make_dicts_formem(nuc,dicts1b,dict_TM,sps,hw)
    Hamil,dictsnt,Chan1b,Chan2bD,Gamma,maxnpq = store_1b2b(sps,dicts1b,dicts,binfo)
    dictTBMEs = dictsnt.dictTBMEs
    MatOp = Matrix{Float64}[]
    if "Rp2" in Operators
        MatOp = [ zeros(Float64,maxnpq,maxnpq) for i=1:2*nthreads()]
    end
    for (i,tnuc) in enumerate(nucs)
        nuc = def_nuc(tnuc,ref,corenuc)
        Z = nuc.Z; N=nuc.N; A=nuc.A
        binfo = basedat(nuc,sntf,hw,emax,ref)
        if i > 1
            recalc_v!(A,dicts)
            update_1b!(binfo,sps,Hamil)
            update_2b!(binfo,sps,Hamil,dictTBMEs,Chan2bD,dicts)
            dictTBMEs = dictsnt.dictTBMEs
        end      
        HFobj = hf_iteration(binfo,HFdata[i],sps,Hamil,dictTBMEs,Chan1b,Chan2bD,Gamma,maxnpq,dWS,to;verbose=verbose,io=io)
        if "Rp2" in Operators
            Op_Rp2 = InitOp(Chan1b,Chan2bD.Chan2b)
            eval_rch_hfmbpt(binfo,Chan1b,Chan2bD,HFobj,Op_Rp2,dWS,MatOp,to;io=io)
        end
    end
    return true
end

function addHCM1b!(Hamil::Operator,HCM::Operator,fac=1.0)
    for pn = 1:2
        Hamil.onebody[pn] += fac .* HCM.onebody[pn]
    end
    return nothing
end

function update_dicts_withHCM!(HCM::Operator,Chan2bD,dicts)
    Chan2b = Chan2bD.Chan2b
    nchan = length(Chan2b)
    @threads for ch = 1:nchan
        tkey = zeros(Int64,4)
        tbc = Chan2b[ch]
        kets = tbc.kets
        J = tbc.J
        Hcm = HCM.twobody[ch]
        pnrank = 2 + div(tbc.Tz,2)
        nkets = length(tbc.kets)
        for i=1:nkets
            bra = kets[i]            
            tbra = @view tkey[1:2] 
            tbra .= bra
            for j=i:nkets
                ket = kets[j]
                tket = @view tkey[3:4] 
                tket .= ket
                tHcm = Hcm[i,j]
                nkey = get_nkey_from_abcdarr(tkey)
                for target in dicts[pnrank][nkey] #[ [totJ,V2b,Vjj,Vjj_2n3n,Vpp*hw,Hcm] ]
                    if J != target[1];continue;end
                    target[6] = tHcm
                end
            end
        end
    end
    return nothing
end

function make_dicts_formem(nuc,dicts1b,dict_TM,sps,hw)
    dicts=[ Dict{Int64,Vector{Vector{Float64}}}() for pnrank=1:3]
    dict_snt2ms = dicts1b.snt2ms    
    A = nuc.A  
    for pnrank = 1:3
        tdict = dict_TM[pnrank]
        target_dict = dicts[pnrank]
        for tkey in keys(tdict)
            val = tdict[tkey]
            a,b,c,d = tkey
            ta = dict_snt2ms[a]; tb = dict_snt2ms[b]          
            tc = dict_snt2ms[c]; td = dict_snt2ms[d]
            oa = sps[ta]; ob = sps[tb]
            oc = sps[tc]; od = sps[td]
            ja = oa.j; jb = ob.j; jc = oc.j; jd = od.j
            key = zeros(Int64,4)
            for tmp in val
                key[1] = ta; key[2] = tb; key[3] = tc; key[4] = td
                totJ,V2bdum,Vjj,Vjj_2n3n,Vpp = tmp                                
                phase_ab = (-1)^(div(ja+jb,2)+totJ+1)
                phase_cd = (-1)^(div(jc+jd,2)+totJ+1)
                flip_ab = ifelse(ta>tb,true,false)
                flip_cd = ifelse(tc>td,true,false)
                phase = 1.0
                if flip_ab; key[1] = tb; key[2] = ta; phase *= phase_ab;end
                if flip_cd; key[3] = td; key[4] = tc; phase *= phase_cd;end
                if key[1] > key[3]
                    k1,k2,k3,k4 = key
                    key[1] = k3; key[2] = k4; key[3] = k1; key[4] = k2
                end
                Vjj *= phase
                Vjj_2n3n *= phase
                Vpp *= phase
                V2b = Vjj + Vjj_2n3n + Vpp*hw/A                
                nkey = get_nkey_from_abcdarr(key)
                t = get(target_dict,nkey,false)
                vcm = 0.0
                if t==false
                    target_dict[nkey] = [[totJ,V2b,Vjj,Vjj_2n3n,Vpp*hw,vcm]]
                else
                    push!(target_dict[nkey],[totJ,V2b,Vjj,Vjj_2n3n,Vpp*hw,vcm])
                end
            end
        end
    end
    return dicts
end

""" 
    prepHFdata(nucs,ref,datatype,corenuc)
Constructor of an array of `hfdata` struct.
"""
function prepHFdata(nucs,ref,datatype,corenuc)
    dnum = length(datatype)    
    Data = hfdata[ ]
    for tnuc in nucs
        nuc = def_nuc(tnuc,ref,corenuc)
        data = [ zeros(Float64,2) for i =1:dnum]
        push!(Data, hfdata(nuc,data,datatype))
    end
    return Data
end 

"""
    def_sps(emax)
Function to define `sps::Vector{SingleParticleState}` from `emax`. 
"""
function def_sps(emax)
    sps = SingleParticleState[ ]
    p_sps = SingleParticleState[ ]
    n_sps = SingleParticleState[ ]
    for temax = 0:emax
        prty = temax % 2
        for l = 0+prty:2:emax
            jmin = 2*l - 1
            jmax = 2*l + 1
            n = div(temax-l,2)
            if n < 0;continue;end
            if jmin < 1;jmin=jmax;end
            for j = jmin:2:jmax
                push!(p_sps,SingleParticleState(n,l,j,-1,[0.0],[false],[false],[false]))
                push!(sps,SingleParticleState(n,l,j,-1,[0.0],[false],[false],[false]))
                push!(n_sps,SingleParticleState(n,l,j,1,[0.0],[false],[false],[false]))
                push!(sps,SingleParticleState(n,l,j,1,[0.0],[false],[false],[false]))
            end
        end
    end
    return sps,p_sps,n_sps
end

"""
    recalc_v!(A,dicts)

Function to re-calculate two-body interaction from snt file for a different mass number.
This is needed because in the readsnt/readsnt_bin function, the interaction part and the kinetic term 
are stored separately to avoid multiple reads of the input file for calculation of multiple nuclei.
"""
function recalc_v!(A,dicts)
    for pnrank = 1:3
        tdict = dicts[pnrank]
        for tkey in keys(tdict)
            tmp = tdict[tkey]            
            for i in eachindex(tmp)
                tmp[i][2] = tmp[i][3] + tmp[i][4] + tmp[i][5]/A + tmp[i][6] *A
            end
        end 
    end
    return nothing
end 

"""
    naive_filling(sps,n_target,emax,for_ref=false)
calculate naive filling configurations by given sps and proton/neutron number (`n_target`)

For some nuclei, especially for heavier ones, carrying out naive filling is ambiguous
(e.g., neutron occupation of 22O can be both 0s1(2),0p1(2),0p3(4),0d5(6) and  0s1(2),0p1(2),0p3(4),1s1(2), 0d3(4)).
In this function, "naive filling" means to try fill orbits with lower ``2n+l`` and then "lower" ``j``.
"""
function naive_filling(sps,n_target,emax,for_ref=false)
    ln = length(sps)
    occ = [ 0.0 for i =1:ln]
    imin = imax = 1
    Nocc = 0
    GreenLight = false
    ofst = 0    
    for e = 0:emax
        j2min = 1
        j2max = 2*e +1
        ncand = sum( [ j2+1 for j2=j2min:2:j2max ])
        cand_j_gt_l = Int64[ ]
        cand_j_lt_l = Int64[ ]
        # GreenLight: whether to fill orbits in a major shell
        if ncand + Nocc <= n_target;GreenLight=true;else;GreenLight=false;end
        for n = ln:-1:1
            if e != 2 * sps[n].n +  sps[n].l;continue;end
            N = sps[n].j +1 
            if GreenLight
                occ[n]=1.0; Nocc += N 
            else
                if sps[n].j == 2*sps[n].l + 1
                    push!(cand_j_gt_l,n)
                else
                    push!(cand_j_lt_l,n)
                end
            end
        end
        if !GreenLight
            fractional = false
            for gt in [true,false]
                cand = ifelse(gt,cand_j_gt_l,cand_j_lt_l)
                imin = ifelse(gt,1,length(cand))
                imax = ifelse(gt,length(cand),1)
                step = ifelse(gt,1,-1)
                for ith = imin:step:imax
                    n = cand[ith]
                    tnocc = sps[n].j + 1
                    if tnocc + Nocc <= n_target
                        occ[n] = 1; Nocc += tnocc
                        if Nocc == n_target; break;end
                    else
                        fractional = true
                        occ[n] = (n_target-Nocc)/tnocc
                        Nocc = n_target
                        break
                    end                
                end
            end
        end 
        if Nocc == n_target && GreenLight 
           break
        end
        ofst += e + 1 
    end
    if Nocc != n_target; println("warn! Nocc");exit();end
    return occ
end

"""
    ini_occ!(pconf,occ_p,nconf,occ_n)

initialize occupation number matrices (```occ_p```&```occ_n```) by naive filling configurations ```pconf```&```nconf```
"""
function ini_occ!(pconf,occ_p,nconf,occ_n)
    for i in eachindex(pconf)
        occ_p[i,i] = pconf[i]
    end
    for i in eachindex(nconf)
        occ_n[i,i] = nconf[i] 
    end    
    return nothing
end

"""
    ReorderHFSPS!(h_p,h_n,Cp,Cn,e1b_p,e1b_n,Chan1b)

"reorder" HF single particle space.
Since we diagonalize the `h_p,h_n` (istead of subblock mat), we need to specify the correspondance between ordering of sps and that of HFSPEs obtained by solving HF eigenvalue problem
"""
function ReorderHFSPS!(h_p,h_n,Cp,Cn,e1b_p,e1b_n,Chan1b; verbose=false)
    for pn = 1:2
        tmp = Chan1b.chs1b[pn]
        tkeys = keys(tmp)
        h = ifelse(pn==1,h_p,h_n)
        evec = ifelse(pn==1,e1b_p,e1b_n)
        C = ifelse(pn==1,Cp,Cn) 
        nonzeros = Dict{Int64,Vector{Int64}}()
        for tkey in tkeys
            idxs = tmp[tkey]
            nidxs = Int64[ ]
            for idx in tmp[tkey]
                nidx = 0
                if pn ==1
                    nidx = div(idx,2) + 1
                else
                    nidx = div(idx,2)
                end
                push!(nidxs,nidx)
            end
            if pn ==1; tkey = div(tkey,2)+1;
            else; tkey = div(tkey,2);end
            for idx in nidxs            
                t = get(nonzeros,tkey,false)
                if t == false
                    nonzeros[tkey] = [idx]
                else
                    if (idx in nonzeros[tkey]) == false
                        push!(nonzeros[tkey],idx)
                    end
                end
                t = get(nonzeros,idx,false)
                if t == false 
                    nonzeros[idx] = [tkey]
                else    
                    if (tkey in nonzeros[idx]) == false
                        push!(nonzeros[idx],tkey)
                    end
                end 
            end
        end
        for tkey = 1:size(C)[1]
            cvec = @views C[:,tkey]
            cvec .= 0.0
            idxs = nonzeros[tkey]
            if length(idxs) == 1
                cvec[tkey] = 1.0
                evec[tkey] = h[tkey,tkey]
            else
                nidxs = sort(idxs)
                sM = @views h[nidxs, nidxs]
                vals,vecs = eigen(sM)
                for (n,idx) in enumerate(nidxs)
                    evec[idx] = vals[n]
                    if vecs[n,n] < 0.0; vecs[:,n] .*= -1.0;end
                    for (m,jdx) in enumerate(nidxs)
                        C[jdx,idx] = vecs[m,n]
                    end
                end
            end
        end
    end
    return nothing
end

function calc_rho!(rho,U,occ,M)   
    BLAS.gemm!('N','T',1.0,occ,U,0.0,M)
    BLAS.gemm!('N','N',1.0,U,M,0.0,rho)
    return nothing
end

"""
    def_holeparticle(Chan1b,occ_p,occ_n,p_sps,n_sps)

define hole/particle space by ```occ_p, occ_n```
"""
function def_holeparticle(Chan1b,occ_p,occ_n,p_sps,n_sps)
    snt2ms = Chan1b.snt2ms
    lp = length(p_sps); ln = length(n_sps)
    particles = [ Int64[ ] , Int64[ ] ]
    holes = [ Int64[ ] , Int64[ ] ]
    for pn = 1:2
        occs = ifelse(pn==1,occ_p,occ_n)
        t_sps = ifelse(pn==1,p_sps,n_sps)
        for i = 1:ifelse(pn==1,lp,ln)
            idx_snt = i + lp*(pn-1)
            msidx = snt2ms[idx_snt]
            t_sps[i].occ[1] = occs[i,i]
            if occs[i,i] == 0.0
                push!(particles[pn],msidx)            
            else
                push!(holes[pn],msidx)
            end
        end
    end
    return holes,particles
end

"""
    calc_Gamma!(Gamma,sps,Cp,Cn,V2,Chan2b,maxnpq)

Function to alculate ``\\Gamma`` (two-body HF interaction).
Note: V3NO from genuine 3NF is supported for ver >= 0.4.0
"""
function calc_Gamma!(Gamma,sps,Cp,Cn,V2,Chan2b,maxnpq,Object_3NF,rho,dWS)
    sps_3b = Object_3NF.sps_3b
    dict_3b_idx = Object_3NF.dict_3b_idx
    v3bme = Object_3NF.v3bme

    l_sps = length(sps)
    E3max = sps_3b.e3max
    nchan = length(Chan2b)
    Ds = [ zeros(Float64,maxnpq,maxnpq) for i =1:nthreads()]
    M  = [ zeros(Float64,maxnpq,maxnpq) for i =1:nthreads()]
    use3NF = Object_3NF.use3BME
    dim_v3 = ifelse(use3NF, maxnpq, 1)
    V3NOs = [ zeros(Float64,dim_v3,dim_v3) for i =1:nthreads()]
    @threads for ch = 1:nchan
        tid = threadid()
        tmp = Chan2b[ch]
        Tz = tmp.Tz; J=tmp.J; kets = tmp.kets
        npq = length(kets)
        D = @view Ds[tid][1:npq,1:npq]; # D .= 0.0
        dim_v3 = ifelse(use3NF, npq, 1)
        V3NO = @view V3NOs[tid][1:dim_v3,1:dim_v3]; V3NO .= 0.0
        v = V2[ch]
        for ib = 1:npq
            i,j = kets[ib]            
            phase_ij = (-1)^( div(sps[i].j+sps[j].j,2) + 1 + J)
            idx_bra1 = div(i,2) + i%2 
            idx_bra2 = div(j,2) + j%2 
            e2bra = (2 * sps[i].n + sps[i].l) + (2 * sps[j].n + sps[j].l)
            for ik = 1:npq
                k,l = kets[ik]
                C1 = ifelse(i%2==1,Cp,Cn)
                C2 = ifelse(j%2==1,Cp,Cn)
                idx_ket1 = div(k,2) + k%2
                idx_ket2 = div(l,2) + l%2
                phase_kl = (-1)^( div(sps[k].j+sps[l].j,2) + 1 + J)
                e2ket = (2 * sps[k].n + sps[k].l) + (2 * sps[l].n + sps[l].l)
                vdtmp = 0.0
                if Tz != 0
                    vdtmp = C1[idx_bra1,idx_ket1] * C2[idx_bra2,idx_ket2]
                    if i!=j
                        vdtmp += C1[idx_bra2,idx_ket1] * C2[idx_bra1,idx_ket2] * phase_ij
                    end
                    if i==j; vdtmp *= sqrt(2.0);end
                    if k==l; vdtmp /= sqrt(2.0);end
                else
                    p_idx_bra = ifelse(i%2==1,idx_bra1,idx_bra2)
                    n_idx_bra = ifelse(i%2==1,idx_bra2,idx_bra1)
                    p_idx_ket = ifelse(k%2==1,idx_ket1,idx_ket2)
                    n_idx_ket = ifelse(k%2==1,idx_ket2,idx_ket1)
                    phase = ifelse(i%2==0,phase_ij,1)           
                    phase *= ifelse(k%2==0,phase_kl,1)
                    vdtmp = Cp[p_idx_bra,p_idx_ket] * Cn[n_idx_bra,n_idx_ket] * phase
                end
                D[ib,ik] = vdtmp
           
                #evaluating_V3NO
                if !use3NF; continue; end
                if ik < ib; continue; end
                v3tmp = 0.0
                for a = 1:l_sps
                    oa = sps[a]
                    ja2 = oa.j
                    ea = 2 * oa.n  + oa.l
                    if ea + e2bra > E3max; continue; end
                    for b = 1:l_sps
                        ob = sps[b]
                        jb2 = ob.j
                        eb = 2 * ob.n + ob.l

                        if oa.l != ob.l || ja2 != jb2 || oa.tz != ob.tz ; continue; end
                        if e2ket + eb > E3max; continue; end

                        J3min = abs(2*J - ja2)
                        J3max = 2*J + ja2
                        for J3 = J3min:2:J3max
                            v3 = get_V3_pn(-1, E3max, v3bme, J, J, J3, i, j, a, k, l, b, sps_3b,dict_3b_idx,dWS)
                            v3tmp += rho[a,b] * (J3+1) *v3
                        end
                    end
                end
                v3tmp /= (2*J + 1)
                if i == j; v3tmp /= sqrt(2.0); end
                if k == l; v3tmp /= sqrt(2.0); end
                V3NO[ib,ik] = v3tmp
                if ib != ik; V3NO[ik,ib] = v3tmp; end
            end
        end
        Gam = Gamma[ch]
        tM  = @views M[threadid()][1:npq,1:npq]
        if use3NF
            BLAS.gemm!('N','N',1.0,v+V3NO,D,0.0,tM)
            BLAS.gemm!('T','N',1.0,D,tM,0.0,Gam)
        else
            BLAS.gemm!('N','N',1.0,v,D,0.0,tM)
            BLAS.gemm!('T','N',1.0,D,tM,0.0,Gam)
        end
       
        # println("ch ",@sprintf("%3i", ch),           
        #         " norm(D) ", @sprintf("%9.3e",norm(D)),
        #         " norm(V2) ", @sprintf("%9.3e",norm(v)),
        #         " norm(V3NO) ", @sprintf("%9.3e",norm(V3NO)),
        #         " norm(Gam) ", @sprintf("%15.7e",norm(Gamma[ch])))
    end
    return nothing
end

function make_symmetric!(mat)
    # this is the brute force way to make array symmetric
    for i=1:size(mat)[1]
        for j=i:size(mat)[1]
            if j!=i; mat[j,i]=mat[i,j];end
        end
    end
    return nothing
end

"""
    add_ch_ket!(ch,iket,tdict) 

add ch & idx for kets in `spaces::space_channel` (pp/hh/etc.)
"""
function add_ch_ket!(ch,iket,tdict) 
    defined = get(tdict,ch,0)
    if defined == 0
        tdict[ch] = [iket]
    else
        push!(tdict[ch],iket)
    end
    return nothing
end

"""
    get_space_chs(sps,Chan2b)

define hole/particle single particle states.
In this function, only the hh/pp/ph (needed for IMSRG) are defined, and other channels will be updated later for target normal ordering or VS-IMSRG flow.
"""
function get_space_chs(sps,Chan2b)    
    hh = Dict{Int64,Vector{Int64}}()
    ph = Dict{Int64,Vector{Int64}}()
    pp = Dict{Int64,Vector{Int64}}()
    cc = Dict{Int64,Vector{Int64}}()
    vc = Dict{Int64,Vector{Int64}}()
    qc = Dict{Int64,Vector{Int64}}()
    vv = Dict{Int64,Vector{Int64}}()
    qv = Dict{Int64,Vector{Int64}}()
    qq = Dict{Int64,Vector{Int64}}()
    for ch in eachindex(Chan2b)
        tbc = Chan2b[ch]
        kets = tbc.kets
        for (ik,ket) in enumerate(kets)
            i,j = ket
            ni = sps[i].occ[1]; nj = sps[j].occ[1]
            if ni + nj == 0.0; 
                add_ch_ket!(ch,ik,pp)
            else
                if ni!=0.0 && nj != 0.0
                    add_ch_ket!(ch,ik,hh)
                else 
                    add_ch_ket!(ch,ik,ph) 
                end
            end
        end       
    end   
    return space_channel(pp,ph,hh,cc,vc,qc,vv,qv,qq)   
end

"""
    getHNO(binfo,tHFdata,E0,p_sps,n_sps,occ_p,occ_n,h_p,h_n,e1b_p,e1b_n,Cp,Cn,V2,Chan1b,Chan2b::tChan2b,Gamma,maxnpq,dict_2b_ch,dWS,to) where{tChan2b <: Vector{chan2b}}

obtain spherical HF solution and calc. MBPT correction (upto 2nd&3rd order) to g.s. energy
"""
function getHNO(binfo,tHFdata,E0,p_sps,n_sps,occ_p,occ_n,h_p,h_n,rho,
                e1b_p,e1b_n,Cp,Cn,V2,Chan1b,Chan2b::Vector{chan2b},Gamma,maxnpq,                
                dict_2b_ch,dWS,Object_3NF::Obj_3BME,to;io=stdout) 
    ## Calc. f (1-body term)
    fp = Cp' * (h_p*Cp); fn = Cn' *(h_n*Cn) # equiv to vals_p/n   
    make_symmetric!(fp); make_symmetric!(fn)

    ## Calc. particle_hole states 
    holes, particles = def_holeparticle(Chan1b,occ_p,occ_n,p_sps,n_sps)
    sps = make_sps_from_pnsps(p_sps,n_sps,Chan1b)
    spaces = get_space_chs(sps,Chan2b)
    modelspace = ModelSpace(p_sps,n_sps,sps,occ_p,occ_n,holes,particles,spaces)
    ## Calc. Gamma (2bchanel matrix element)    
    @timeit to "Gamma" calc_Gamma!(Gamma,sps,Cp,Cn,V2,Chan2b,maxnpq,Object_3NF,rho,dWS)
    EMP2 = HF_MBPT2(binfo,modelspace,fp,fn,e1b_p,e1b_n,Chan2b,Gamma;io=io)
    EMP3 = HF_MBPT3(binfo,modelspace,e1b_p,e1b_n,Chan2b,dict_2b_ch,dWS,Gamma,to;io=io)
    exists = get(ame2020data,binfo.nuc.cnuc,false)   
    Eexp = 0.0
    if exists==false
        println(io,"E_HF ", @sprintf("%12.5f",E0), 
        "  E_MBPT(3) = ",@sprintf("%12.4f",E0+EMP2+EMP3),"  Eexp: Not Available")
    else
        Eexp = - binfo.nuc.A * ame2020data[binfo.nuc.cnuc].BE/1000.0
        println(io,"E_HF ", @sprintf("%12.5f",E0),
        "  E_MBPT(3) = ",@sprintf("%12.4f",E0+EMP2+EMP3),"  Eexp: "*@sprintf("%12.3f", Eexp))  
    end
    println("")
    tmp = tHFdata.data
    E = tmp[1]
    E[1] = E0+EMP2+EMP3; E[2] = Eexp
    H0 = Operator([E0],[fp,fn],Gamma,[true],[false])
    return HamiltonianNormalOrdered(H0,E0,EMP2,EMP3,Cp,Cn,e1b_p,e1b_n,modelspace)
end

function get_rho!(rho, rho_p,rho_n)
    @assert size(rho_p) == size(rho_n) "rho_p and rho_n must have the same size"
    Dim = size(rho)[1]
    rho[1:2:Dim,1:2:Dim] .= rho_p
    rho[2:2:Dim,2:2:Dim] .= rho_n
    return nothing
end

"""
    hf_iteration(binfo,tHFdata,sps,Hamil,dictTBMEs,Chan1b,Chan2bD,Gamma,maxnpq,dWS,to;itnum=100,verbose=false,HFtol=1.e-14,inttype="snt")

solve HF equation

This function returns object with HamiltonianNormalOrdered (HNO) struct type, which contains...
- `E0,EMP2,EMP3` HF energy and its MBPT corrections
- `fp/fn::Matrix{Float64}` one-body int.
- `Gamma:: Vector{Matrix{Float64}}` two-body int.
"""
function hf_iteration(binfo,tHFdata,sps,Hamil,dictTBMEs,Chan1b,Chan2bD,Gamma,maxnpq,dWS,to;
                      Object_3NF = "", itnum=300,verbose=false,HFtol=1.e-14,io=stdout,E0cm=0.0)                      
    Chan2b = Chan2bD.Chan2b; dict_2b_ch = Chan2bD.dict_ch_JPT
    dim1b = div(length(sps),2)
    mat1b = zeros(Float64,dim1b,dim1b)
    p1b = Hamil.onebody[1]
    n1b = Hamil.onebody[2]
    V2 = Hamil.twobody
    nuc = binfo.nuc; emax=binfo.emax
    Z = nuc.Z; N = nuc.N   
    abcd = zeros(Int64,4)
    p_sps,n_sps = get_pn_sps(sps)
    occ_p = zeros(Float64,dim1b,dim1b); occ_n = zeros(Float64,dim1b,dim1b)
    EHFs = [ zeros(Float64,6) for i=1:2]
    pconf = naive_filling(p_sps,Z,emax);nconf = naive_filling(n_sps,N,emax)
    ini_occ!(pconf,occ_p,nconf,occ_n)
    rho_p = copy(mat1b); Cp = copy(mat1b); Up = copy(mat1b);for i=1:dim1b;Up[i,i]=occ_p[i,i];end
    rho_n = copy(mat1b); Cn = copy(mat1b); Un = copy(mat1b);for i=1:dim1b;Un[i,i]=occ_n[i,i];end
    calc_rho!(rho_p,Up,occ_p,Cp);calc_rho!(rho_n,Un,occ_n,Cn)
    e1b_p = zeros(Float64,dim1b); e1b_n = zeros(Float64,dim1b)
    Vt_pp = copy(mat1b); Vt_nn = copy(mat1b); Vt_pn = copy(mat1b); Vt_np = copy(mat1b)
    calc_Vtilde(sps,Vt_pp,Vt_nn,Vt_pn,Vt_np,rho_p,rho_n,dictTBMEs,abcd,Chan1b)
    V3tilde = zeros(Float64,0,0)
    rho = zeros(Float64, dim1b*2, dim1b*2)
    @timeit to "eval_V3NO" if Object_3NF.use3BME        
        get_rho!(rho, rho_p,rho_n)
        V3tilde = zeros(Float64, length(sps), length(sps))
        eval_V3NO!(sps,V3tilde,rho,rho_p,rho_n,Object_3NF)        
    end
    h_p = copy(mat1b); h_n = copy(mat1b)
    update_FockMat!(h_p,p1b,p_sps,h_n,n1b,n_sps,Vt_pp,Vt_nn,Vt_pn,Vt_np,Object_3NF,V3tilde)
    calc_Energy(rho_p,rho_n,p1b,n1b,p_sps,n_sps,Vt_pp,Vt_nn,Vt_pn,Vt_np,EHFs,V3tilde,Object_3NF;verbose=true)

    for it = 1:itnum        
        ## diagonalize proton/neutron 1b hamiltonian
        valsp,vecsp = eigen(h_p); valsn,vecsn = eigen(h_n)

        ## Update 1b density matrix
        Up .= vecsp; Un .= vecsn
        ReorderHFSPS!(h_p,h_n,Up,Un,valsp,valsn,Chan1b)
        calc_rho!(rho_p,Up,occ_p,Cp);calc_rho!(rho_n,Un,occ_n,Cn)     
   
        ## Re-evaluate tilde(V) and Fock matrix
        calc_Vtilde(sps,Vt_pp,Vt_nn,Vt_pn,Vt_np,rho_p,rho_n,dictTBMEs,abcd,Chan1b)
        eval_V3NO!(sps,V3tilde,rho,rho_p,rho_n,Object_3NF)        
        
        update_FockMat!(h_p,p1b,p_sps,h_n,n1b,n_sps,Vt_pp,Vt_nn,Vt_pn,Vt_np,Object_3NF,V3tilde)
        calc_Energy(rho_p,rho_n,p1b,n1b,p_sps,n_sps,Vt_pp,Vt_nn,Vt_pn,Vt_np,EHFs,V3tilde,Object_3NF;verbose=verbose)
        if HF_conv_check(EHFs;tol=HFtol)
            #print("HF converged @ $it  \t")
            printEHF(EHFs[1])
            valsp,vecsp = eigen(h_p); valsn,vecsn = eigen(h_n)
            e1b_p .= valsp;e1b_n .= valsn; Up .= vecsp; Un .= vecsn    
            ReorderHFSPS!(h_p,h_n,Up,Un,e1b_p,e1b_n,Chan1b)       
            calc_rho!(rho_p,Up,occ_p,Cp);calc_rho!(rho_n,Un,occ_n,Cn)    
            get_rho!(rho, rho_p, rho_n)

            Cp .= vecsp; Cn .= vecsn
            ReorderHFSPS!(h_p,h_n,Cp,Cn,e1b_p,e1b_n,Chan1b)
            break
        end
        tnorm = norm(Up'*Up-Matrix{Float64}(I, dim1b,dim1b),Inf)
        if tnorm > 1.e-10;println("Unitarity check: res. norm(p) $tnorm");end
    end
    ## HNO: get normal-ordered Hamiltonian
    E0 = EHFs[1][1] - E0cm
    @timeit to "getHNO" HFobj = getHNO(binfo,tHFdata,E0,p_sps,n_sps,occ_p,occ_n,h_p,h_n,rho,e1b_p,e1b_n,Cp,Cn,V2,Chan1b,Chan2b,Gamma,maxnpq,dict_2b_ch,dWS,Object_3NF,to;io=io)
    return HFobj
end

function eval_V3NO!(sps,V3tilde,rho,rho_p,rho_n,Object_3NF)
    get_rho!(rho, rho_p, rho_n)
    dict_idx_me3j_to_snt = Object_3NF.dict_idx_to_snt
    V3tilde .= 0.0
    v3monopole = Object_3NF.v3monopole
    v = [ zeros(Float64,size(V3tilde)[1],size(V3tilde)[2]) for i = 1:nthreads()]    
    keylist = collect(keys(v3monopole))
    @threads for idx in eachindex(keylist)
        tkey = keylist[idx]
        a,c,i,b,d,j = unhash_key6j(tkey)
        rho_ab = rho[a,b]; rho_cd = rho[c,d]
        
        a = dict_idx_me3j_to_snt[a]
        b = dict_idx_me3j_to_snt[b]
        c = dict_idx_me3j_to_snt[c]
        d = dict_idx_me3j_to_snt[d]
        i = dict_idx_me3j_to_snt[i]
        j = dict_idx_me3j_to_snt[j]

        v3tmp = v[threadid()]
        v3tmp[i,j] += rho_ab * rho_cd * v3monopole[tkey]
    end
    for i = 1:nthreads()
        V3tilde .+= v[i]
    end
    V3tilde .+= transpose(V3tilde) - Diagonal(V3tilde)
    
    return nothing
end

"""
    printEHF(Es)

print HF energy and its break down ```Es=[E1b,E2bpp,E2bnn,E2bpn,E3b]```
"""
function printEHF(Es)
    println("E: ", @sprintf("%12.6f", Es[1]),"  = E1b ", @sprintf("%9.5f", Es[2]),"  + E2b ", @sprintf("%9.5f", Es[3]+Es[4]+Es[5]),
            "   ( "*@sprintf("%9.3f", Es[3])*@sprintf("%9.3f", Es[4]),@sprintf("%9.3f", Es[5]),"), + E3b ", @sprintf("%9.5f", Es[6]))
end

function calc_Energy(rho_p,rho_n,p1b,n1b,p_sps,n_sps,Vt_pp,Vt_nn,Vt_pn,Vt_np,Es,V3tilde,Object_3NF;verbose=false)
    ## 1-body part
    lp = size(p1b)[1];ln = size(n1b)[1]
    ep_1b = 0.0; en_1b = 0.0
    for alph =1:lp
        Na = p_sps[alph].j * 1.0 + 1.0        
        for beta = 1:lp
            ep_1b += p1b[alph,beta] * rho_p[alph,beta] * Na
        end        
    end
    for alph =1:ln
        Na = n_sps[alph].j * 1.0 + 1.0
        for beta = 1:ln
            en_1b += n1b[alph,beta] * rho_n[alph,beta] *Na
        end
    end
    E1b = ep_1b + en_1b

    ## 2-body part
    E2b = 0.0; E2bpp = 0.0; E2bpn = 0.0; E2bnn = 0.0   
    for i = 1:lp
        for j=1:lp
            if rho_p[i,j] == 0.0;continue;end
            E2bpp += 0.5 * rho_p[i,j] *Vt_pp[i,j]
            E2bpn += 0.5 * rho_p[i,j] *Vt_pn[i,j]
        end        
    end
    for i = 1:ln
        for j=1:ln
            if rho_n[i,j] == 0.0;continue;end
            E2bnn += 0.5* rho_n[i,j] *Vt_nn[i,j]
            E2bpn += 0.5* rho_n[i,j] *Vt_np[i,j]
        end        
    end
    E2b = E2bpp + E2bpn + E2bnn

    ## 3-body (NO2B) part
    E3b = 0.0
    if size(V3tilde)[1] != 0
        dict_to_me3j = Object_3NF.dict_idx_to_me3j
        rho = zeros(Float64,2*lp,2*lp)
        get_rho!(rho, rho_p,rho_n)
        for i = 1:size(V3tilde)[1]
            idx_i = div(i,2) + ifelse(i%2==1,1,0)
            Jfactor = p_sps[idx_i].j + 1
            for j = 1:size(V3tilde)[2]
                rho_ij = rho[i,j]
                idx_i = idx_j = -1
                if i%2 == j%2 == 1
                    idx_i = div(i,2) + 1
                    idx_j = div(j,2) + 1
                else
                    idx_i = div(i,2) 
                    idx_j = div(j,2) 
                end
                V3ij = V3tilde[dict_to_me3j[i],dict_to_me3j[j]]

                if V3ij  != 0.0
                    E3b += 1/6 * rho_ij * V3ij  * Jfactor
                    if rho_ij != 0.0
                        V2 = 0.0
                        if i%2 == j%2 == 1
                            idx_i = div(i,2) + 1
                            idx_j = div(j,2) + 1
                            N = p_sps[idx_j].j + 1
                            V2 += (Vt_pn[idx_i,idx_j] + Vt_pp[idx_i,idx_j]) / N 
                        else
                            idx_i = div(i,2) 
                            idx_j = div(j,2) 
                            N = p_sps[idx_j].j + 1
                            V2 += (Vt_np[idx_i,idx_j] + Vt_nn[idx_i,idx_j] ) / N 
                        end                        
                    end
                end
            end
        end
    end

    E = E1b + E2b + E3b
    if verbose
        println("E:", @sprintf("%15.6f", E),"E1b:", @sprintf("%15.6f", E1b),
                "E2b:", @sprintf("%15.6f", E2b)," E3b ", @sprintf("%15.6f", E3b))
    end
    Es[2] .= Es[1]
    Es[1] .= [ E,E1b,E2bpp,E2bpn,E2bnn, E3b]
    return nothing
end

function HF_conv_check(EHFs;tol=1.e-9)
    old = EHFs[2]; new = EHFs[1]
    if (abs(old[1] - new[1]) < tol) 
        return true
    else
        return false
    end
end

"""
$(SIGNATURES)

Functon updating Fock matrix.
Since the ``F_{ij}``
"""
function update_FockMat!(h_p,p1b,p_sps,h_n,n1b,n_sps,Vt_pp,Vt_nn,Vt_pn,Vt_np, Object_3NF,V3tilde)
    lp = size(h_p)[1]; ln = size(h_n)[1]
    h_p .= p1b; h_n .= n1b
    ## for proton
    for i = 1:lp
        Ni = p_sps[i].j + 1.0
        for j = 1:lp # Vpp
            h_p[i,j] += (Vt_pp[i,j]+Vt_pn[i,j]) / Ni
        end
    end
    ## for neutron
    for i = 1:ln
        Ni = n_sps[i].j+ 1.0
        for j = 1:ln # Vnn
            h_n[i,j] += (Vt_nn[i,j]+Vt_np[i,j]) / Ni
        end
    end

    dict_idx_me3j_to_snt = Object_3NF.dict_idx_to_snt
    if Object_3NF.use3BME
        @assert lp == ln "3BME is only implemented for lp=ln model space"
        Dim = 2 * lp
        for i = 1:2:Dim # i (proton)
            idx_i = div(dict_idx_me3j_to_snt[i],2) + 1
            Ni = p_sps[idx_i].j + 1.0
            for j = 1:2:Dim # j (proton)
                idx_j = div(dict_idx_me3j_to_snt[j],2) + 1
                v3tmp = V3tilde[i,j]
                h_p[idx_i,idx_j] += V3tilde[i,j] / 2#Ni
            end
            for j = 2:2:Dim # j (neutron)
                idx_j = div(dict_idx_me3j_to_snt[j],2)
                h_p[idx_i,idx_j] += V3tilde[i,j] / 2#Ni
            end
        end
        for i = 2:2:Dim # i (neutron)
            idx_i = div(dict_idx_me3j_to_snt[i],2)
            Ni = n_sps[idx_i].j + 1.0
            for j = 1:2:Dim # j = proton)
                idx_j = div(dict_idx_me3j_to_snt[j],2)+1
                h_n[idx_i,idx_j] += V3tilde[i,j] / 2#Ni
            end
            for j = 2:2:Dim # j (neutron)
                idx_j = div(dict_idx_me3j_to_snt[j],2)
                h_n[idx_i,idx_j] += V3tilde[i,j] / 2#Ni
            end
        end
    end
    return nothing
end

function calc_Vtilde(sps,Vt_pp,Vt_nn,Vt_pn,Vt_np,rho_p,rho_n,dictTBMEs,tkey,Chan1b;verbose=false)
    symfac = 1.0/3
    dim1b = size(Vt_pp)[1]
    dict_pp = dictTBMEs[1];dict_pn = dictTBMEs[2];dict_nn = dictTBMEs[3]    
    Vt_pp .= 0.0; Vt_nn .= 0.0; Vt_pn .= 0.0;Vt_np .= 0.0
    Chan1b_p,Chan1b_n = Chan1b.chs1b
    for idx_i = 1:dim1b
        i = 2*(idx_i-1) + 1
        for idx_j = idx_i:dim1b
            j = 2*(idx_j-1)+1
            if !(j in Chan1b_p[i]);continue;end
            ji = sps[i].j; ji = sps[j].j; 
            ## tilde(Vp) from Vpp
            for idx_a = 1:dim1b
                a = 2*(idx_a-1) + 1
                for idx_b = idx_a:dim1b
                    b = 2*(idx_b-1) + 1
                    if !(b in Chan1b_p[a]);continue;end
                    rho_ab = rho_p[idx_a,idx_b]
                    tkey[1] = i; tkey[3] = j;tkey[2] = a; tkey[4] = b
                    if a < i; tkey[2] = i; tkey[4] = j;tkey[1] = a; tkey[3] = b; end
                    vmono,vmono2n3n = dict_pp[tkey]
                    Vt_pp[idx_i,idx_j] += rho_ab * (vmono + vmono2n3n * symfac)
                    # if rho_ab != 0.0; println("pp: i $i j $j  a $a b $b rho $rho_ab key $tkey vmono ",vmono/(sps[i].j+1));end
                    if a!=b
                        if a < i                       
                            tkey[3] = a; tkey[1] = b
                        else
                            tkey[4] = a; tkey[2] = b
                        end
                        vmono,vmono2n3n = dict_pp[tkey]
                        Vt_pp[idx_i,idx_j] += rho_ab * (vmono + vmono2n3n * symfac)
                        # if rho_ab !=0.0; println("pp: i $i j $j  a $a b $b rho $rho_ab key $tkey vmono ",vmono/(sps[i].j+1));end
                    end
                end
            end
            Vt_pp[idx_j,idx_i] = Vt_pp[idx_i,idx_j]
            ## tilde(Vp) from Vpn
            for idx_a= 1:dim1b
                a = 2*idx_a
                for idx_b = idx_a:dim1b
                    b = 2*idx_b
                    if !(b in Chan1b_n[a]);continue;end
                    rho_ab = rho_n[idx_a,idx_b] 
                    tkey[1] = i;tkey[3] = j; tkey[2] = a; tkey[4] = b
                    vmono,vmono2n3n = dict_pn[tkey]
                    Vt_pn[idx_i,idx_j] += rho_ab * (vmono + vmono2n3n * symfac) 
                    #if rho_ab != 0.0; println("pn: i $i j $j  a $a b $b rho $rho_ab key $tkey vmono ",vmono/(sps[i].j+1));end
                    if a!=b
                        tkey[1] = i; tkey[3] = j; tkey[2] = b; tkey[4] = a
                        vmono,vmono2n3n = dict_pn[tkey]
                        Vt_pn[idx_i,idx_j] += rho_ab *  (vmono + vmono2n3n* symfac)
                        #if rho_ab != 0.0; println("pn: i $i j $j  a $a b $b rho $rho_ab key $tkey vmono ",vmono/(sps[i].j+1));end
                    end
                end
            end
            Vt_pn[idx_j,idx_i] = Vt_pn[idx_i,idx_j]
        end
    end
    for idx_i = 1:dim1b
        i = 2*idx_i 
        for idx_j = idx_i:dim1b
            j = 2*idx_j
            if !(j in Chan1b_n[i]);continue;end
            ## tilde(Vn) from Vnn
            for idx_a = 1:dim1b
                a = 2*idx_a 
                for idx_b = idx_a:dim1b
                    b = 2*idx_b 
                    if !(b in Chan1b_n[a]);continue;end
                    rho_ab = rho_n[idx_a,idx_b]
                    tkey[1] = i; tkey[3] = j;tkey[2] = a; tkey[4] = b                    
                    if a < i; tkey[2] = i; tkey[4] = j;tkey[1] = a; tkey[3] = b; end
                    vmono,vmono2n3n = dict_nn[tkey]
                    Vt_nn[idx_i,idx_j] += rho_ab * (vmono + vmono2n3n* symfac) 
                    if a!=b
                        if a < i                       
                            tkey[3] = a; tkey[1] = b
                        else
                            tkey[4] = a; tkey[2] = b
                        end
                        vmono,vmono2n3n = dict_nn[tkey]
                        Vt_nn[idx_i,idx_j] += rho_ab * (vmono + vmono2n3n* symfac)
                    end
                end
            end
            Vt_nn[idx_j,idx_i] = Vt_nn[idx_i,idx_j]
            ## tilde(Vn) from Vnp
            for idx_a= 1:dim1b
                a = 2*(idx_a-1) + 1
                for idx_b = idx_a:dim1b
                    b = 2*(idx_b-1) + 1
                    if !(b in Chan1b_p[a]);continue;end
                    rho_ab = rho_p[idx_a,idx_b] 
                    tkey[1] = a ; tkey[3] = b
                    tkey[2] = i; tkey[4] = j
                    vmono,vmono2n3n = dict_pn[tkey]
                    Vt_np[idx_i,idx_j] += rho_ab * (vmono + vmono2n3n*symfac)
                    if a!=b
                        tkey[1] = b; tkey[3] = a
                        tkey[2] = i; tkey[4] = j
                        vmono,vmono2n3n = dict_pn[tkey]
                        Vt_np[idx_i,idx_j] += rho_ab * (vmono + vmono2n3n*symfac)
                    end
                end
            end
            Vt_np[idx_j,idx_i] = Vt_np[idx_i,idx_j]
        end
    end
    return nothing
end
