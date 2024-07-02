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

function call_SVD(X, r, method)
    if method == "Arpack"
        Z = svds(X; nsv=r)[1]
        return Z.U, Z.S, Z.Vt
    elseif method == "KrylovKit"
        m, n = size(X)
        vals, lvecs, rvevs, info  = svdsolve(X, m, r)
        U_r = zeros(Float64, m, r)
        S_r = zeros(Float64, r)
        Vt_r = zeros(Float64, r, n)
        S_r .= vals[1:r]
        for i = 1:r
            U_r[:,i] .= lvecs[i]
            Vt_r[i,:] .= rvevs[i]
        end
        return U_r, S_r, Vt_r
    else
        if method != "full"
            println("$method is not supported! svd in LinearAlgebra is called instead.")
        end
        Z = svd(X)
        if r <= rank(X) 
            return Z.U[:,1:r], Z.S[1:r], Z.Vt[1:r,:]
        else
            return Z.U, Z.S, Z.Vt            
        end
    end
end

function get_DMD_operator(U_r, S_r, Vt_r, Y)
    S_r_inv = diagm( 1.0 ./ S_r)
    Z = BLAS.gemm('T', 'N', 1.0, Vt_r, S_r_inv)
    YZ = BLAS.gemm('N','N',1.0, Y, Z)
    Atilde = BLAS.gemm('T', 'N', 1.0, U_r, YZ)
    return Atilde
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
    tnorm = max(norm(Y - Yapprox,Inf), norm(Y - Yapprox,2))
    return tnorm
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
function extrapolate_DMD(x_start, U_r, Atilde, s_pred, fn_exact, s_end, ds, nuc, inttype, emax, oupdir; verbose=false)

Function to perform the DMD extrapolation.
The time evolution (IMSRG flow) of original data `x` is approximated by the time evolution of 
`z` in the latent space, and then the approximated data `x'` is written to a HDF5 file.
"""
function extrapolate_DMD(x_start, U_r, Atilde, s_pred, fn_exact, s_end, ds, nuc, inttype, emax, oupdir; verbose=false)
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
            while s < s_target
                z_k .= z_new
                BLAS.gemv!('N', 1.0, Atilde, z_k, 0.0, z_new)
                s += ds
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
        if verbose
            check_stationarity(z1_r, z_k, z_new, Atilde)
        end
    end
    return nothing
end

function write_dmdvec_hdf5(vec_in,s,E_imsrg,nuc,inttype,emax,oupdir)
    vec = zeros(Float64,length(vec_in)+2)
    vec[1] = s
    vec[2] = E_imsrg
    vec[3:end] .= vec_in
    fname = oupdir*"omega_dmdvec"*ifelse(inttype=="","","_"*inttype)*"_e$(emax)_$(nuc)_s"*strip(@sprintf("%6.2f",s))*".h5"
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

function util_SVD(X, Y, r_max, tol_svd, method, allow_fullSVD, to; dont_care_stationarity=false)
    @assert method == "full" || method == "Arpack" || method == "KrylovKit" "method must be full, KrylovKit, or Arpack"
    fullrank = rank(X)
    r_used = min(r_max,fullrank)
    println("fullrank $fullrank r_used $r_used")
    if allow_fullSVD || method == "full"
        @timeit to "full SVD" U,Sigma,Vt = call_SVD(X, r_used, "full")
        r_used = min(r_max, fullrank, sum(Sigma .> tol_svd))
        print_vec("checking full Sigma[1:$r_used]...", Sigma[1:r_used];ine=true)
    end
    @timeit to "SVD($method)" U_r, S_r, Vt_r = call_SVD(X,r_used,method)

    redright = true
    u = @view U_r[:,1:r_used]
    s = @view S_r[1:r_used]
    vt = @view Vt_r[1:r_used,:]
    Atilde = get_DMD_operator(u,s,vt,Y)
    while redright && !dont_care_stationarity 
        Atilde = get_DMD_operator(u,s,vt,Y)
        evals_Atilde = eigvals(Atilde)
        redright = !check_stationarity_Atilde(evals_Atilde)
        print_vec("checking Sigma...", S_r[1:r_used];ine=true)
        if redright == false; break; end
        println("truncation at rank <= $r_used does not give converged results")
        r_used -= 1
        @assert r_used > 0 "Aborted because the snapshots don't give converged resutls. Please modify smin/smax of snapshots."
        u = @view U_r[:,1:r_used]
        s = @view S_r[1:r_used]
        vt = @view Vt_r[1:r_used,:]
    end
    print_vec("singular values r <= $r_used ", S_r[1:r_used];ine=true)
    return fullrank, r_used, u, Atilde
end

function check_packages_SVD(X)
    r = min(10,rank(X))
    m, n = size(X)
    SVD = svd(X)
    s = SVD.S    
    s_K = svdsolve(X, m, r)[1]
    s_A = svds(X; nsv=r)[1].S

    @assert issorted(s,rev=true) "s must be sorted in descending order"

    println("Comparing the singular values by different packages....")
    print_vec("full SVD (LinearAlgebra)", s[1:r];ine=true)
    print_vec("trun.SVD (KrylovKit)    ", s_K[1:r];ine=true)
    print_vec("trun.SVD (Arpack)       ", s_A[1:r];ine=true)
    return nothing
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
- `tol_svd::Float64`: tolerance for the singular values for the truncated SVD, too small value gives noisy behavior
- `inttype::String`: interaction type, which is used for the filename of the output data
- `methodSVD::String`: method for the truncated SVD.
    - `"full"`: full SVD (LinearAlgebra)
    - `"Arpack"`: truncated SVD svds Arpack.jl
    - `"KrylovKit"`: truncated SVD using `svdsolve[1]` in KrylovKit.jl, this may be not appropriate for the current problem, so do not use it here.
- `oupdir::String`: output directory for the emulated dmdvec data
- `is_show::Bool`: if `true`, the TimerOutput result is shown
- `debugmode::Bool`: if `true`, the packages for the SVD are checked
- `dont_care_stationarity::Bool`: if `true`, the stationarity of the Atilde is not checked
"""
function dmd_main(emax, nuc, fns, r_max, smin, ds;s_pred=Float64[],fn_exact=String[],
                  allow_fullSVD=true,tol_svd=1e-6,inttype="",
                  methodSVD="Arpack", oupdir="flowOmega/",is_show=true, debugmode=false,
                  dont_care_stationarity=true)
    to = TimerOutput()
    if !isdir("flowOmega")
        println("dir. flowOmega is created!")
        mkdir("flowOmega")
    end
    println("Trying to perform DMD....")

    #  construct snapshot matrices X and Y
    s_end, X,Y = constructor_snapshot_matrix(fns)
    println("Snapshot from smin $smin s_end $s_end ds $ds")

    if debugmode
        check_packages_SVD(X)
    end
    
    fullrank, r, U_r, Atilde = util_SVD(X, Y, r_max, tol_svd, methodSVD, allow_fullSVD, to; dont_care_stationarity=dont_care_stationarity) 
    evals_Atilde = eigvals(Atilde)
    check_stationarity_Atilde(evals_Atilde)
    plot_Atilde_eigvals(evals_Atilde, emax, nuc, smin, s_end, ds, inttype, fullrank, r)

    if is_show
        show(to); println()
    end
    # extrapolation
    extrapolate_DMD(Y[:,end],U_r, Atilde, s_pred, fn_exact, s_end, ds, nuc, inttype, emax, oupdir)

    return nothing
end

function check_stationarity_Atilde(evals; tol=1.e-2)
    tf = all(abs.(evals) .<= 1.0+tol)
    println("eigenvalues of Atilde ", evals)
    println("stationary?  may be... ", tf)
    return tf
end

function plot_Atilde_eigvals(evals, emax, nuc, smin, s_end, ds, inttype, fullrank, r)
    if !isdir("pic") 
        mkdir("pic")
    end
    fn = "pic/Atilde_eigvals_e$(emax)_$(nuc)_smin$(smin)_send$(s_end)_ds$(ds)_rank$(r)_of_$(fullrank).pdf"
    # make figure square
    p = plot(size=(300,300))

    # draw unit circle
    θ = LinRange(0,2π,100)
    plot!(cos.(θ),sin.(θ),lw=2,alpha=0.8,label=:none)

    # plot eigenvalues .<= 1.0
    inon = evals[findall(abs.(evals) .<= 1.0)]
    plot!(real(inon),imag(inon),st=:scatter,label=L"|z| \leq 1", alpha=0.8)
    
    out = evals[findall(abs.(evals) .> 1.0)]
    if length(out) > 0
        plot!(real(out),imag(out),st=:scatter,label=L"|z| > 1",alpha=0.8)
    end

    xlabel!(L"\mathrm{Re} (z)")
    ylabel!(L"\mathrm{Im} (z)")
    savefig(p,fn)
    
    return nothing
end