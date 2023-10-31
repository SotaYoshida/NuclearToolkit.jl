
"""
read_fvec_hdf5(fname)

Function to read a flattened vector from a HDF5 file.
Note that the first and second element of the vector are the value of `s` and `Es`, respectively,
and the rest of the vector is the flattened form of the Operator object.
"""
function read_fvec_hdf5(fname)
    io = h5open(fname,"r")
    s = read(io,"s")
    Es = read(io,"Es")
    dim = read(io,"dim")
    vec = zeros(Float64,dim+2)
    vec[1] = s
    vec[2] = Es
    idx_s = 3
    vec_p1b = read(io,"p1b"); vec[idx_s:idx_s+length(vec_p1b)-1] .= vec_p1b; idx_s += length(vec_p1b) 
    vec_n1b = read(io,"n1b"); vec[idx_s:idx_s+length(vec_n1b)-1] .= vec_n1b; idx_s += length(vec_n1b) 
    vec_pp2b = read(io,"pp2b"); vec[idx_s:idx_s+length(vec_pp2b)-1] .= vec_pp2b; idx_s += length(vec_pp2b)
    vec_pn2b = read(io,"pn2b"); vec[idx_s:idx_s+length(vec_pn2b)-1] .= vec_pn2b; idx_s += length(vec_pn2b)
    vec_nn2b = read(io,"nn2b"); vec[idx_s:idx_s+length(vec_nn2b)-1] .= vec_nn2b; idx_s += length(vec_nn2b)    
    close(io)
    return vec  
end

"""
write_fvec_hdf5(binfo,fvec,dict_if_idx_for_hdf5,s,Es;label="omega")

Function to write a flattened vector to a HDF5 file.
Note that the first and second element of the vector are the value of `s` and `Es`, respectively,
and the rest of the vector is the flattened form of the Operator object.
The correspondance between the index of the flattened vector and the Operator object is given by `dict_if_idx_for_hdf5`.
"""
function write_fvec_hdf5(binfo,fvec,dict_if_idx_for_hdf5,s,Es;label="omega")
    pid = getpid()
    fname = "flowOmega/$(label)_vec_$pid"*binfo.nuc.cnuc*"_s"*strip(@sprintf("%6.2f",s))*".h5"
    io = h5open(fname,"w")
    dim = length(fvec) 
    idx_i,idx_j = dict_if_idx_for_hdf5["p1b"]; vec_p1b = fvec[idx_i:idx_j]
    idx_i,idx_j = dict_if_idx_for_hdf5["n1b"]; vec_n1b = fvec[idx_i:idx_j]
    idx_i,idx_j = dict_if_idx_for_hdf5["pp2b"]; vec_pp2b = fvec[idx_i:idx_j]
    idx_i,idx_j = dict_if_idx_for_hdf5["pn2b"]; vec_pn2b = fvec[idx_i:idx_j]
    idx_i,idx_j = dict_if_idx_for_hdf5["nn2b"]; vec_nn2b = fvec[idx_i:idx_j]
    write(io,"s",s)
    write(io,"Es",Es)
    write(io,"dim",dim-2)
    write(io,"p1b",vec_p1b)
    write(io,"n1b",vec_n1b)
    write(io,"pp2b",vec_pp2b)
    write(io,"pn2b",vec_pn2b)
    write(io,"nn2b",vec_nn2b)
    close(io)
    return nothing
end

"""
make_Op_from_flattenvector(fvec,similar_to::Operator,dict_idx_flatvec_to_op,ovwrite_zerobody=false)
    
Function to convert a flattened vector to an Operator object.
Note that the zerobody part is zero unless the `ovwrite_zerobody` is set to `true`.
"""
function make_Op_from_flattenvector(fvec,similar_to::Operator,dict_idx_flatvec_to_op,ovwrite_zerobody=false)
    Op = deepcopy(similar_to); aOp!(Op,0.0)
    if ovwrite_zerobody
        Op.zerobody[1] = fvec[2]
    end
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

"""
update_Op_with_fvec!(fvec,Op,dict_idx_flatvec_to_op)    
    
Function to convert a flattened vector to an Operator object.
Note that the zerobody part is zero unless the `ovwrite_zerobody` is set to `true`.
"""
function update_Op_with_fvec!(fvec,Op,dict_idx_flatvec_to_op,ovwrite_zerobody=false)    
    if ovwrite_zerobody
        Op.zerobody[1] = fvec[2]
    end
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

"""
get_fvec_from_Op(s, Op::Operator,dict_idx_op_to_flatvec, dict_idx_flatvec_to_op)
    
It returns a flattened vector from an Operator object.
Note that the first and second element of the vector are the value of `s` and `Op.zerobody[1]`, respectively.
"""
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

"""
get_fvec_from_Op!(s, vec,Op::Operator,dict_idx_op_to_flatvec, dict_idx_flatvec_to_op)

The destructive version of `get_fvec_from_Op`.
It overwrites the flatten vector by the elements of the Operator object.
"""
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
function get_non0omega_idxs(HFobj::HamiltonianNormalOrdered,Omega::Operator,Chan2b,mode="normal",verbose=false)
    @assert (mode == "normal" || mode=="valence") "unsupported mode for get_non0omega_idxs: $mode"
    dim1b = size(Omega.onebody[1])[1]
    p_sps = HFobj.modelspace.p_sps
    n_sps = HFobj.modelspace.n_sps

    dict_idx_op_to_flatvec=Dict{Vector{Int64},Int64}()
    dict_idx_flatvec_to_op=Dict{Int64,Vector{Int64}}()
    dict_if_idx_for_hdf5 = Dict{String,Vector{Int64}}()

    hit = hitnon0 = 2 # since we keep s(E0) explicitly as the first (second) element of flatten vector
    idx_s = hit + 1
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
                            #println("1bch-idx pn $pn i $i $oi j $j $oj")
                        end
                    else
                        # to be extended to VS-IMSRG
                    end
                end
            end
        end
        target = ifelse(pn==1,"p1b","n1b")
        dict_if_idx_for_hdf5[target] = [idx_s,hitnon0]
        idx_s = hitnon0 + 1
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

    Tz = -4
    ppidx = pnidx = nnidx = 0
    for ch = 1:nch
        tbc = Chan2b[ch]
        @assert tbc.Tz >= Tz "$(tbc.Tz) > $Tz should never happen"
        Tz = tbc.Tz
        dim = div(tbc.nkets*(tbc.nkets+1),2)
        if Tz == -2;ppidx += dim; end
        if Tz ==  0;pnidx += dim; end
        if Tz ==  2;nnidx += dim; end
    end
    dict_if_idx_for_hdf5["pp2b"] = [idx_s, idx_s+ppidx-1]; idx_s += ppidx
    dict_if_idx_for_hdf5["pn2b"] = [idx_s, idx_s+pnidx-1]; idx_s += pnidx
    dict_if_idx_for_hdf5["nn2b"] = [idx_s, idx_s+nnidx-1]; idx_s += nnidx
    if verbose
        println("dim. flatvec $(hit) # of nonzero element (other than E(s)) $hitnon0")
    end
    return dict_idx_op_to_flatvec, dict_idx_flatvec_to_op, dict_if_idx_for_hdf5
end


"""
svd_Op_twobody(s,Op::Operator,Chan2b;verbose=true,max_rank=20)

Function to perform SVD of two-body operator `Op`.
"""
function svd_Op_twobody(s,Op::Operator,Chan2b;verbose=true,max_rank=20)
    println("SVD @s=$s")
    dimfull = dimtr = 0
    for ch in eachindex(Chan2b)
        tbc = Chan2b[ch]
        J = tbc.J; P =tbc.prty; Tz = tbc.Tz
        dim = tbc.nkets
        if dim == 0; continue;end
        mat = Op.twobody[ch]
        fullrank = rank(mat)
        SVD = LinearAlgebra.svd(mat)
        if verbose 
            print("ch ",@sprintf("%5i",ch), " JPT=($J,$P,$Tz)\t fullrank ", @sprintf("%5i",fullrank), "/ dim= ",@sprintf("%5i",dim),"  ")
            #print_vec("singular values", SVD.S)
            tol = 1.e-8
            hit = 0
            dimfull += div(dim*(dim+1),2)
            for trank = 1:fullrank
                U = SVD.U; Sig = deepcopy(SVD.S); Vt = SVD.Vt
                Sig[trank+1:end] .= 0.0
                SV = Diagonal(Sig)* Vt
                Vtilde = BLAS.gemm('N','N',1.0,U,SV)
                tnorm = norm(mat-Vtilde,2)
                if tnorm < tol
                    hit += 1
                    println("rank=",@sprintf("%4i",trank), " norm(V-V') ",@sprintf("%12.5e",norm(mat-Vtilde,2)))
                    dimtr += 2*trank*dim + trank
                    # if trank != fullrank 
                    #     trank = fullrank 
                    #     U = SVD.U; Sig = deepcopy(SVD.S); Vt = SVD.Vt
                    #     Sig[trank+1:end] .= 0.0
                    #     SV = Diagonal(Sig)* Vt
                    #     Vtilde = BLAS.gemm('N','N',1.0,U,SV)
                    #     tnorm = norm(mat-Vtilde,2)
                    #     println("rank=",@sprintf("%4i",trank), " norm(V-V') ",@sprintf("%12.5e",norm(mat-Vtilde,2)))
                    # end
                    break
                end
            end
            if hit == 0; println("");end
        end        
    end
    println("dimfull $dimfull dimtr $dimtr")
    return nothing
end