function ann_imsrg()

end

function read_fvec_bin(fname)
    io = open(fname,"r")
    s = read(io,Float64)
    Es = read(io,Float64)
    dim = read(io,Int64)
    vec = zeros(Float64,dim+2)
    vec[1] = s
    vec[2] = Es
    vec[3:end] .= [read(io,Float64) for i = 1:dim]
    close(io)
    return vec
end

function write_fvec_bin(binfo,fvec,istep,s,Es;label="omega")
    pid = getpid()
    nw = istep
    fname = "flowOmega/$(label)_vec_$pid"*binfo.nuc.cnuc*"_$nw.bin"
    io = open(fname,"w")
    dim = length(fvec) 
    write(io,Float64(s))
    write(io,Float64(Es))
    write(io,Int64(dim-2))
    write(io,fvec[3:end])
    close(io)
    return nothing
end

function make_Op_from_flattenvector(fvec,similar_to::Operator,dict_idx_flatvec_to_op)
    Op = deepcopy(similar_to); aOp!(Op,0.0)
    #Op.zerobody[1] = fvec[2]
    for tkey in keys(dict_idx_flatvec_to_op)
        ch, i, j = dict_idx_flatvec_to_op[tkey]
        target = Op.onebody[1]
        if ch == 0
            target = Op.onebody[2]
        elseif ch  > 0 
            target = Op.twobody[ch]
        end
        target[i,j] = fvec[tkey]
        target[j,i] =-fvec[tkey]
    end
    return Op
end 

function update_Op_with_fvec!(fvec,Op,dict_idx_flatvec_to_op)    
    #Op.zerobody[1] = fvec[2]
    aOp!(Op,0.0)
    for tkey in keys(dict_idx_flatvec_to_op)
        ch, i, j = dict_idx_flatvec_to_op[tkey]
        target = Op.onebody[1]
        if ch == 0
            target = Op.onebody[2]
        elseif ch > 0 
            target = Op.twobody[ch]
        end
        target[i,j] =  fvec[tkey]
        target[j,i] = -fvec[tkey]
    end
    return nothing
end

function get_fvec_from_Op(s, Op::Operator,dict_idx_op_to_flatvec, dict_idx_flatvec_to_op)
    vec = zeros(Float64,2 + length(keys(dict_idx_flatvec_to_op)))
    vec[1] = s
    vec[2] = Op.zerobody[1]
    for tkey in keys(dict_idx_op_to_flatvec)
        fidx = dict_idx_op_to_flatvec[tkey]
        ch, i, j = tkey
        target = Op.onebody[1]
        if ch == 0
            target = Op.onebody[2]
        elseif ch > 0
            target = Op.twobody[ch]
        end
        vec[fidx] = target[i,j] 
    end
    return vec
end

function get_fvec_from_Op!(s, vec,Op::Operator,dict_idx_op_to_flatvec, dict_idx_flatvec_to_op)
    vec .*= 0.0
    vec[1] = s
    vec[2] = Op.zerobody[1]
    for tkey in keys(dict_idx_op_to_flatvec)
        fidx = dict_idx_op_to_flatvec[tkey]
        ch, i, j = tkey
        target = Op.onebody[1]
        if ch == 0
            target = Op.onebody[2]
        elseif ch > 0
            target = Op.twobody[ch]
        end
        vec[fidx] = target[i,j] 
    end
    return nothing
end

"""
Function to make Operator flatten vector for ANN calculations, which are to be performed with flattened vectors.
mode: can be "normal" or "valence"

For VS-IMSRG, there is a possibility to work to deal with only `vv` space.
"""
function get_non0omega_idxs(HFobj::HamiltonianNormalOrdered,Omega::Operator,mode="normal",verbose=false)
    @assert (mode == "normal" || mode=="valence") "unsupported mode for get_non0omega_idxs: $mode"
    dim1b = size(Omega.onebody[1])[1]
    p_sps = HFobj.modelspace.p_sps
    n_sps = HFobj.modelspace.n_sps

    dict_idx_op_to_flatvec=Dict{Vector{Int64},Int64}()
    dict_idx_flatvec_to_op=Dict{Int64,Vector{Int64}}()
    hit = hitnon0 = 2 # since we keep s(E0) explicitly as the first (second) element of flatten vector
    for pn = 1:2
        ch = pn-2
        ms = ifelse(pn==1,p_sps,n_sps)
        for i = 1:dim1b
            oi = ms[i]
            for j=1:dim1b
                oj = ms[j]
                if i < j 
                    hit += 1
                    if mode == "normal"                    
                        if oi.l == oj.l && oi.j == oj.j && oi.tz == oj.tz
                            hitnon0 += 1
                            dict_idx_op_to_flatvec[[ch,i,j]] = hitnon0
                            dict_idx_flatvec_to_op[hitnon0] = [ch,i,j]
                        end
                    else
                        # to be extended to VS-IMSRG
                    end
                end
            end
        end
    end    
    nch = length(Omega.twobody)
    for ch = 1:nch
        Mat2b = Omega.twobody[ch]
        nket = size(Mat2b)[1]
        for i = 1:nket
            for j = 1:nket
                if i <= j
                    hit += 1
                    #     ppidx = get(HFobj.modelspace.spaces.pp,ch,Int64[])
                    #     hhidx = get(HFobj.modelspace.spaces.hh,ch,Int64[])                    
                    hitnon0 += 1
                    dict_idx_op_to_flatvec[[ch,i,j]] = hitnon0
                    dict_idx_flatvec_to_op[hitnon0] = [ch,i,j]
                end
            end
        end
    end

    if verbose
        println("dim. flatvec $(hit) # of nonzero element (other than E(s)) $hitnon0")
    end
    return dict_idx_op_to_flatvec, dict_idx_flatvec_to_op
end

function imsrg_flow_check(nucs,sntf,hw,emax_calc;verbose=false,Operators=String["Rp2"],doIMSRG=false,valencespace=[],corenuc="",ref="nucl",return_obj=false,oupfn="",debugmode=0,fn_params="optional_parameters.jl",Hsample=false,emulator=true)
    @assert isfile(sntf) "sntf:$sntf is not found!"
    to = TimerOutput()
    io = select_io(false,"",nucs;use_stdout=true,fn=oupfn)
    chiEFTparams = init_chiEFTparams(;io=nothing)
    HFdata = prepHFdata(nucs,ref,["E"],corenuc)
    @timeit to "prep dWS2n" dWS = prep_dWS2n(chiEFTparams,to;emax_calc=emax_calc)
    @timeit to "read" begin        
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

        MatOp = [ zeros(Float64,maxnpq,maxnpq) for i=1:2*nthreads()]
    end
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
            HFobj = hf_iteration(binfo,HFdata[i],sps,Hamil,dictsnt.dictTBMEs,Chan1b,Chan2bD,Gamma,maxnpq,dWS,to;verbose=verbose,io=io,E0cm=E0cm) 
        end
        IMSRGobj = imsrg_main(binfo,Chan1b,Chan2bD,HFobj,dictsnt,dWS,valencespace,Operators,MatOp,to;fn_params=fn_params,debugmode=debugmode,Hsample=Hsample,emulator=emulator)
        Aold = A

    end    
    return true
end