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
                    break
                end
            end
            if hit == 0; println("");end
        end        
    end
    println("dimfull $dimfull dimtr $dimtr")
    return nothing
end

function constructor_snapshot_matrix(fns)
    num_snapshots = length(fns)
    fvec_dim = length(read_fvec_hdf5(fns[1])) -2 

    X = zeros(Float64,fvec_dim,num_snapshots-1)
    Y = zeros(Float64,fvec_dim,num_snapshots-1)
    for (i,fname) in enumerate(fns[1:end-1])
        fvec = read_fvec_hdf5(fname)[3:end]
        X[:,i] .= fvec
        if i > 1
            Y[:,i-1] .= fvec
        end
    end
    Y[:,end] .= read_fvec_hdf5(fns[end])[3:end]
    s_end = split(split(fns[end],"_s")[end],".h5")[1]
    return parse(Float64,s_end),X,Y
end

function get_DMD_operator(X,Y,r)
    Z = svds(X; nsv=r)[1]
    U_r, S_r, Vt_r = Z.U, Z.S, Z.Vt

    S_r_inv = diagm( 1.0 ./ S_r)
    Z = BLAS.gemm('T', 'N', 1.0, Vt_r, S_r_inv)
    YZ = BLAS.gemm('N','N',1.0, Y, Z)
    Atilde = BLAS.gemm('T', 'N', 1.0, U_r, YZ)

    return U_r, Atilde
end

function check_DMD_norm(X, Y, r, U_r, Atilde; verbose=false)
    Y_latent = zeros(Float64, r, size(Y)[2])
    x1 = X[:,1]
    x1_r = BLAS.gemv('T', 1.0, U_r, x1)
    x_k = zeros(Float64,r)
    x_new = zeros(Float64,r) .+ x1_r
    for k = 1:size(Y)[2]
        x_k .= x_new
        BLAS.gemv!('N', 1.0, Atilde, x_k, 0.0, x_new)
        Y_latent[:,k] .= x_new
        
    end
    Yapprox = BLAS.gemm('N', 'N', 1.0, U_r, Y_latent)
    if verbose 
        for k = 1:size(Y)[2]
            println("k $k normvec ", norm(Yapprox[:,k]-Y[:,k]))
            print_vec("Yap", Yapprox[1:10,k])
            print_vec("Y  ", Y[1:10,k])
        end
    end
    println("norm(Y-Yapprox,Inf) = ", @sprintf("%10.4e",norm(Y - Yapprox,Inf)), " Fro. ", @sprintf("%10.4e",norm(Y - Yapprox,2)))
end

function check_stationarity(z, z_k, z_pred,  Atilde, itnum=1000)
    tf = true
    z .= z_pred
    z_k .= z_pred
    norms = zeros(Float64,itnum)
    for it = 1:itnum        
        BLAS.gemv!('N', 1.0, Atilde, z, 0.0, z_k)
        norms[it] = norm(z_k - z_pred)
    end
    println("norms $norms")
    return tf
end

"""
function extrapolate_DMD(x_start, U_r, Atilde, s_pred, fn_exact, s_end, ds, nuc, inttype, emax, oupdir)

Function to perform the DMD extrapolation.
The time evolution (IMSRG flow) of original data `x` is approximated by the time evolution of 
`z` in the latent space, and then the approximated data `x'` is written to a HDF5 file.
"""
function extrapolate_DMD(x_start, U_r, Atilde, s_pred, fn_exact, s_end, ds, nuc, inttype, emax, oupdir)
    @assert length(s_pred) == length(fn_exact) "s_pred and fn_exact must have the same length"
    if length(s_pred) > 0
        r = size(U_r)[2]
        z1_r = BLAS.gemv('T', 1.0, U_r, x_start)
        z_k = zeros(Float64,r)
        z_new = zeros(Float64,r) .+ z1_r
        for ith = 1:length(s_pred)
            s_target = s_pred[ith]
            z_k .= 0.0
            z_new .= z1_r
            s = s_end
            while true
                z_k .= z_new
                BLAS.gemv!('N', 1.0, Atilde, z_k, 0.0, z_new)
                s += ds
                if s >= s_target; break; end
            end
            x_pred = BLAS.gemv('N', 1.0, U_r, z_k)
            E_imsrg = 0.0
            if isfile(fn_exact[ith])
                fvec_inf = read_fvec_hdf5(fn_exact[ith])
                s_file = fvec_inf[1]
                @assert s == s_file "s $s must be equal to s_file $s_file for $(fn_exact[ith])"
                E_imsrg = fvec_inf[2]
                x_inf = fvec_inf[3:end]
                println("s =  ",@sprintf("%6.2f", s),"   ||x'-x||  ", @sprintf("%10.4e", norm(x_inf-x_pred)), "   ", @sprintf("%10.4e",norm(x_inf-x_pred,Inf)))
            end
            write_dmdvec_hdf5(x_pred,s,E_imsrg,nuc,inttype,emax,oupdir)
        end
        check_stationarity(z1_r, z_k, z_new, Atilde)
    end
    return nothing
end

function write_dmdvec_hdf5(vec_in,s,E_imsrg,nuc,inttype,emax,oupdir)
    vec = zeros(Float64,length(vec_in)+2)
    vec[1] = s
    vec[2] = E_imsrg
    vec[3:end] .= vec_in
    fname = oupdir*"omega_dmdvec_$(inttype)_e$(emax)_$(nuc)_s"*strip(@sprintf("%6.2f",s))*".h5"
    io = h5open(fname,"w")
    write(io,"vec",vec)
    close(io)
    return nothing
end

function read_dmdvec_hdf5(fn)
    io = h5open(fn,"r")
    vec = read(io,"vec")
    close(io)
    return vec
end

"""
main API for DMD

# Arguments
- `emax::Int64`: maximum energy for the IMSRG calculation, which is used only for the filename of the emulated fvec data
- `nuc::String`: nucleus name
- `fns::Vector{String}`: filenames of the snapshot data
- `trank::Int64`: specified largest rank of truncated SVD
- `smin::Float64`: starting value of `s` for the training data
- `ds::Float64`: step size of `s` for the training data

# Optional arguments
- `s_pred::Vector{Float64}`: values of `s` for the extrapolation
- `fn_exact::Vector{String}`: filenames of the exact data for the extrapolation, which must have the same length as `s_pred`
- `allow_fullSVD::Bool`: if `true`, the full SVD is performed and the rank is determined by `trank` and the tolerance `tol_svd`
- `tol_svd::Float64`: tolerance for the singular values for the truncated SVD
- `inttype::String`: interaction type, which is used for the filename of the output data
"""
function dmd_main(emax, nuc, fns, trank, smin, ds;s_pred=Float64[],fn_exact=String[],
                  allow_fullSVD=true,tol_svd=1e-7,inttype="",oupdir="flowOmega/")
    if !isdir("flowOmega")
        println("dir. flowOmega is created!")
        mkdir("flowOmega")
    end
    println("Trying to perform DMD....")

    #  construct snapshot matrices X and Y
    s_end, X,Y = constructor_snapshot_matrix(fns)
    println("Snapshot from smin $smin s_end $s_end ds $ds")

    # truncated SVD using Arpack.jl and construct tilde(A)
    fullrank = rank(X)
    r = trank
    if allow_fullSVD
        SVD = svd(X)
        sigma_full = SVD.S
        r = min(r, fullrank, sum(sigma_full .> tol_svd))
        println("fullrank $(fullrank) rank $r")
        print_vec("singular values", sigma_full[1:r];ine=true)
    end
    U_r, Atilde = get_DMD_operator(X,Y,r)    
    check_DMD_norm(X, Y, r, U_r, Atilde)

    # extrapolation
    extrapolate_DMD(Y[:,end],U_r, Atilde, s_pred, fn_exact, s_end, ds, nuc, inttype, emax, oupdir)

    return nothing
end
