struct cfps
    chs::Vector{Vector{Matrix{Int64}}}
    cfpdims::Vector{Vector{Vector{Int64}}}
    cfpvals::Vector{Vector{Matrix{Float64}}}
end

"""
# key structure for integ_D [Ebra,Eket]=>[e1,e1p,e2,e2p]=>[K,K1,X]
"""
struct inner_integ_nonlocal
    integ_E::Dict{Int64,Float64}
    integ_D::Dict{Int64, Dict{Int64, Dict{Int64, Dict{Int64, Dict{Int64, Float64}}}}}
end
"""
``\\int d\\pi_1 d\\pi_2 ``
"""
function inner_nonlocal(params,to)
    N3max = params.N3max
    Rnls_r,Rnls_p = prep_Rnls(params)
    gKXs = prep_gKX(params)
    ps = params.meshpoints.ps
    ws = params.meshpoints.wps
    inner_E = Dict{Int64,Float64}()
    for E = 0:N3max
        for n1 = 0:div(E,2)
            n2 = div(E - 2*n1,2)
            e1 = 2*n1 
            e2 = 2*n2
            if E != e1+e2; continue; end
            tkey1 = get_nkey3(e1,n1,0)
            tkey2 = get_nkey3(e2,n2,0)
            R1s = Rnls_p[tkey1]
            R2s = Rnls_p[tkey2]
            ret = 0.0
            for (imesh1,w1) in enumerate(ws)
                p1 = ps[imesh1]
                R1 = R1s[imesh1]
                for (imesh2,w2) in enumerate(ws)
                    p2 = ps[imesh2]
                    R2 = R2s[imesh2]
                    freg = nonlocal_regulator_3NF(p1,p2,2,params.Lambda3NF;type="Navratil")
                    tsum = w1 * w2 * p1 * p2 * R1 * R2 * (-1)^(n1+n2) * freg
                    ret += tsum
                end
            end
            tkey = get_nkey3(E,n1,n2)
            inner_E[tkey] = ret            
        end
    end      
    inner_D = Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}}}()
    for E_bra = 0:N3max
        for E_ket = E_bra:N3max
            Ekey = get_nkey2(E_bra,E_ket)
            inner_D[Ekey] = Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}}()
            target = inner_D[Ekey]
            for e1 = 0:2:E_bra
                e2 = E_bra - e1 
                n = div(e1,2)
                l = 0
                tkey1 = get_nkey3(e1,n,l); Rn0s = Rnls_p[tkey1]    
                for N = 0:div(e2,2)
                    L = e2 - 2*N
                    tkey2 = get_nkey3(e2,N,L); RNLs = Rnls_p[tkey2]
                    for e1p = 0:2:E_ket
                        e2p = E_ket - e1p                         
                        e4key = get_nkey4(e1,e2,e1p,e2p)
                        tmp = get(target,e4key,nothing)
                        if tmp == nothing
                            target[e4key] = Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}()
                        end
                        tdict = target[e4key]
                        np = div(e1p,2)
                        lp = 0
                        tkey1p = get_nkey3(e1p,np,lp); Rn0s_p = Rnls_p[tkey1p]
                        for Np = 0:div(e2p,2)
                            Lp = e2p - 2*Np
                            tkey2p = get_nkey3(e2p,Np,Lp); RNLs_p = Rnls_p[tkey2p]
                            key_bra = get_nkey4(n,l,N,L)
                            key_ket = get_nkey4(np,lp,Np,Lp)
                            tmp = get(tdict,key_bra,nothing)
                            if tmp == nothing;tdict[key_bra] = Dict{Int64,Dict{Int64,Float64}}();end                       
                            if get(tdict[key_bra],key_ket,nothing)==nothing;tdict[key_bra][key_ket] = Dict{Int64,Float64}();end
                                
                            ret = tdict[key_bra][key_ket]
                            for K = 0:2:2
                                for K1 = 0:K
                                    Xmin = max(abs(L-K1),abs(Lp-(K-K1)))
                                    Xmax = min(L+K1,Lp+K-K1)
                                    for X = Xmin:Xmax
                                        tkey = get_nkey3(K,K1,X)
                                        ret[tkey] = 0.0
                                    end
                                end
                            end
                            for K = 0:2:2
                                for K1 = 0:K
                                    Xmin = max(abs(L-K1),abs(Lp-(K-K1)))
                                    Xmax = min(L+K1,Lp+K-K1)
                                    for X = Xmin:Xmax
                                        gKXdict = gKXs[get_nkey2(K,X)]
                                        KK1Xkey = get_nkey3(K,K1,X)
                                        for (i_pi1,w1) in enumerate(ws)
                                            pi1 = ps[i_pi1]
                                            Rn0 = Rn0s[i_pi1] * pi1
                                            for (i_pi2,w2) in enumerate(ws)
                                                pi2 = ps[i_pi2]
                                                RNL = RNLs[i_pi2] * pi2
                                                Fbra = nonlocal_regulator_3NF(pi1,pi2,2,params.Lambda3NF;type="Navratil")
                                                for (i_pi1p,w1p) in enumerate(ws)
                                                    pi1p = ps[i_pi1p]
                                                    Rn0p = Rn0s_p[i_pi1] * pi1p
                                                    for (i_pi2p,w2p) in enumerate(ws)
                                                        pi2p = ps[i_pi2p]
                                                        RNLp = RNLs_p[i_pi2p] * pi2
                                                        Fket = nonlocal_regulator_3NF(pi1p,pi2p,2,params.Lambda3NF;type="Navratil")
                                                        pi2key = get_nkey2(i_pi2,i_pi2p)
                                                        integ = Rn0 * RNL * Fbra * Rn0p * RNLp *Fket * pi2^K1 * pi2p^(K-K1) * gKXdict[pi2key]
                                                        integ *= w1 * w2 * w1p * w2p
                                                        ret[KK1Xkey] += integ 
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return inner_integ_nonlocal(inner_E,inner_D)
end

function gKX(K,X,p,pp,us,ws)
    Mpi_fm = mean(mpis) / hc

        ret = 0.0
    for (i,u) in enumerate(us)
        w = ws[i]
        px = Legendre(X,u)
        pc2 =  (p^2+pp^2-2*p*pp*u)
        nume = pc2^(1-K/2)
        deno = 2/3 * pc2 + Mpi_fm^2
        ret += w * px * nume / deno
    end
    return 0.5 * ret
end

function prep_gKX(params)
    Nmax = params.N3max
    ps = params.meshpoints.ps
    us = params.meshpoints.us
    ws = params.meshpoints.wus
    gKXs = Dict{Int64,Dict{Int64,Float64}}()
    for K = 0:2:2
        for X = 0:Nmax+2
            KXkey = get_nkey2(K,X) 
            gKXs[KXkey] = Dict{Int64,Float64}()
            for (ip,p) in enumerate(ps)
                for (ipp,pp) in enumerate(ps)
                    tkey = get_nkey2(ip,ipp)
                    gKXs[KXkey][tkey] = gKX(K,X,p,pp,us,ws)
                end
            end            
        end
    end
    return gKXs
end

struct inner_integ_local
    integ_c1::Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}}}
    integ_c34::Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}}}
    integ_D::Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}}}
    integ_E::Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}}}
end

function check_construct_integdict(mydict,Ekey,e4key,key_bra=-1,key_ket=-1;level=1)
    @assert 1 <= level <= 2 "Option Error!: level must be 1,2!!"
    if level ==1
        tmp = get(mydict[Ekey],e4key,nothing)
        if tmp==nothing; mydict[Ekey][e4key] = Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}();end
    elseif level == 2
        @assert key_bra != -1
        tmp = get(mydict[Ekey][e4key],key_bra,nothing)
        if tmp == nothing
            mydict[Ekey][e4key][key_bra] = Dict{Int64,Dict{Int64,Float64}}()
        end
        @assert key_ket != -1
        tmp = get(mydict[Ekey][e4key][key_bra],key_ket,nothing)
        if tmp == nothing
            mydict[Ekey][e4key][key_bra][key_ket] = Dict{Int64,Float64}()
        end
    end
    return nothing
end

"""
function to calculate integrals in local Jacobi HO three-body matrix elements
"""
function inner_local(params,to)
    N3max = params.N3max
    Rnls_r,Rnls_p = prep_Rnls(params)
    ps = params.meshpoints.ps
    ws = params.meshpoints.wps     
    @timeit to "Zs" Z0s,Z0Xs,fKs,fKXs = prep_Z0s_Z0Xs_fKs(params)
    Xmax_tpe = params.j3max + 1 + 2 + 2
    #println("Xmax_tpe $Xmax_tpe N3max $N3max L3max $(params.L3max) L3max*2 + 2(K2) + 2(K3) $(params.L3max*2+2+2) ")
  
    integ_E = Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}}}()
    integ_D = Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}}}()
    integ_c1 = Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}}}()
    integ_c34 = Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}}}()
  
    @threads for n_thre = 1:4
        target = integ_E
        if n_thre == 2; target = integ_D;end
        if n_thre == 3; target = integ_c1;end
        if n_thre == 4; target = integ_c34;end

        for E_bra = 0:N3max
            for E_ket = E_bra:N3max
                Ekey = get_nkey2(E_bra,E_ket)
                target[Ekey] = Dict{Int64,Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}}()
                for e1 = 0:E_bra
                    e2 = E_bra - e1 
                    for n1 = 0:div(e1,2)
                        l1 = e1 - 2*n1
                        tkey1 = get_nkey3(e1,n1,l1)
                        R1s = Rnls_r[tkey1]
                        for n2 = 0:div(e2,2)
                            l2 = e2 - 2*n2
                            tkey2 = get_nkey3(e2,n2,l2)
                            R2s = Rnls_r[tkey2]
                            for e1p = 0:E_ket
                                e2p = E_ket - e1p 
                                e4key = get_nkey4(e1,e2,e1p,e2p)
                                check_construct_integdict(target,Ekey,e4key)
                                for n1p = 0:div(e1p,2)
                                    l1p = e1p - 2*n1p
                                    tkey1p = get_nkey3(e1p,n1p,l1p)
                                    R1ps = Rnls_r[tkey1p]
                                    for n2p = 0:div(e2p,2)
                                        l2p = e2p - 2*n2p
                                        tkey2p = get_nkey3(e2p,n2p,l2p)
                                        R2ps = Rnls_r[tkey2p]
                                        key_bra = get_nkey4(n1,l1,n2,l2)
                                        key_ket = get_nkey4(n1p,l1p,n2p,l2p)
                                        integ_inner_nnlo(n_thre,target,Ekey,e4key,key_bra,key_ket,l1,l1p,l2,l2p,Xmax_tpe,
                                                         ps,ws,R1s,R1ps,R2s,R2ps,Z0s,Z0Xs,fKs,fKXs)
                                    end
                                end
                            end
                        end
                    end 
                end
            end
        end
    end
    return inner_integ_local(integ_c1,integ_c34,integ_D,integ_E)
end

function integ_inner_nnlo(n_thre,target,Ekey,e4key,key_bra,key_ket,l1,l1p,l2,l2p,Xmax_tpe,
                         ps,ws,R1s,R1ps,R2s,R2ps,Z0s,Z0Xs,fKs,fKXs)
    check_construct_integdict(target,Ekey,e4key,key_bra,key_ket;level=2)
    ret = target[Ekey][e4key][key_bra][key_ket]
    if n_thre == 1  # Contact
        Xmin = max(abs(l1-l1p),abs(l2-l2p))
        Xmax = min(l1+l1p,l2+l2p)
        for X = Xmin:Xmax  
            tsum = 0.0
            for (imesh1,w1) in enumerate(ws)
                R1 = R1s[imesh1]; R1p = R1ps[imesh1]
                Z0 = Z0s[imesh1]
                for (imesh2,w2) in enumerate(ws)
                    R2  = R2s[imesh2]; R2p = R2ps[imesh2]
                    Z0X = Z0Xs[get_nkey3(X,imesh1,imesh2)]
                    tsum += w1 * w2 * R1 * R1p * R2 * R2p * Z0 * Z0X
                end
            end
            ret[X] = tsum
        end
    elseif n_thre == 2 # OPE
        Xmin_D = abs(l2-l2p)
        Xmax_D = l2 + l2p
        for K = 0:2:2
            for X = Xmin_D:Xmax_D
                tsum = 0.0
                for (imesh1,w1) in enumerate(ws)
                    R1 = R1s[imesh1]; R1p = R1ps[imesh1]
                    fK = fKs[get_nkey2(K,imesh1)]
                    for (imesh2,w2) in enumerate(ws)
                        R2  = R2s[imesh2]; R2p = R2ps[imesh2]
                        Z0X = Z0Xs[get_nkey3(X,imesh1,imesh2)]
                        tsum +=  w1 * w2 * R1 * R1p * R2 * R2p * fK * Z0X
                    end
                end                                        
                ret[get_nkey2(K,X)] = tsum                                        
            end
        end
    elseif n_thre == 3 # TPE: c1 term
        for K3 = 0:1
            for X = 0:Xmax_tpe
                tsum = 0.0
                for (imesh1,w1) in enumerate(ws)
                    R1 = R1s[imesh1]; R1p = R1ps[imesh1];p1 = ps[imesh1]
                    f1  = fKs[get_nkey2(1,imesh1)]
                    for (imesh2,w2) in enumerate(ws)
                        R2  = R2s[imesh2]; R2p = R2ps[imesh2];p2 = ps[imesh2]
                        f1X = fKXs[get_nkey4(1,X,imesh1,imesh2)]
                        tsum += w1 * w2 * R1 * R1p * R2 * R2p * f1 * f1X * (p1/sqrt(2.0))^K3 * (p2*sqrt(3/2))^(1-K3)
                    end
                end
                ret[get_nkey2(K3,X)] = tsum
            end
        end
    elseif n_thre == 4 # TPE: c3&c4 term
        for X = 0:Xmax_tpe
            for K1 = 0:2:2
                for K2 = 0:2:2                
                    for K3 = 0:K2                                                
                        tsum = 0.0
                        for (imesh1,w1) in enumerate(ws)
                            R1 = R1s[imesh1]; R1p = R1ps[imesh1];p1 = ps[imesh1]
                            fK1 = fKs[get_nkey2(K1,imesh1)]
                            for (imesh2,w2) in enumerate(ws)
                                R2  = R2s[imesh2]; R2p = R2ps[imesh2];p2 = ps[imesh2]
                                fK2X = fKXs[get_nkey4(K2,X,imesh1,imesh2)]
                                tsum += w1 * w2 * R1 * R1p * R2 * R2p * fK1 * fK2X * (p1/sqrt(2.0))^K3 * (p2*sqrt(3/2))^(K2-K3)
                            end
                        end
                        ret[get_nkey4(K1,K2,K3,X)] = tsum
                    end
                end
            end
        end        
    end
    return nothing
end

function prep_inner_integ(params3N,to)
    regulator = params3N.regulator  
    tfunc = ifelse(regulator=="local",inner_local,inner_nonlocal)
    return tfunc(params3N,to)
end

function Calc_3NF_in_JacobiHO_coorrdinate(params3N,LECs,to)
    print("Calculating inner integrals... ")
    @timeit to "prep_inner" inner_integ = prep_inner_integ(params3N,to)
    println("=> Done!")
    labframe = true
    if labframe
        #calculation of 3NF in laboratory frame(???)
        print("Calculating CFPs")
        @timeit to "chJPT" JPTs,cfps = calc_channel_JPT(params3N,to)
        println("=> Done!")
        print("Calculating JacobiHO matrix elements")
        @timeit to "calc JacobiHO_3NF" JacobiHO_3NF(LECs,params3N,JPTs,cfps,inner_integ)
        println("=> Done!")
        @timeit to "read JacobiHO_3NF" read_JacobiHO_3NF(params3N,JPTs)
    else
        #solve three-body system
    end    
    return nothing
end

"""
calculate three-body channels {J,P,T}
"""
function calc_channel_JPT(params,to;verbose=false)
    if !isdir("./A3files"); mkdir("A3files"); end
    j3max = params.j3max
    JPTs = Vector{Int64}[ ]
    n = 0
    for t = 1:2:3
        for j = 1:2:j3max
            for parity = 1:-2:-1
                n += 1
                push!(JPTs, [j,parity,t])
            end 
        end
    end    
    nch = (j3max +1) * 2
    @assert n== nch "warn! n $n must be nch $nch in function calc_channel_JPT"
    cfpchs = Vector{Matrix{Int64}}[ ] 
    cfpdims = [ Vector{Int64}[] for ich = 1:nch]
    cfpvals = [ Matrix{Float64}[] for ich = 1:nch]
    for ich = 1:nch
        cfpdim_ch = cfpdims[ich]
        cfpval_ch = cfpvals[ich]
        JPT = JPTs[ich]      
        j,p,t = JPT 
        if verbose 
            println("ich ",@sprintf("%4i",ich), "   J =",@sprintf("%3i", j),
                    "   P =",@sprintf("%3i", p),"   T =",@sprintf("%3i", t))
        end
        nlsjts = get_dim_cfp!(cfpdim_ch,cfpval_ch,JPT,params,to)
        push!(cfpchs,nlsjts)
        if !check_binexists(JPT,params,"cfp")
            set_cfps(cfpdim_ch,cfpval_ch,JPT,nlsjts,params,to)
        else
            read_cfp_bin!(cfpdim_ch,cfpval_ch,JPT,params,to)
        end
    end 
    return JPTs,cfps(cfpchs,cfpdims,cfpvals)
end

function make_binname(JPT,params,bintype;dirname="./A3files/")
    j,p,t=JPT
    Nmax = params.N3max
    fname = dirname*"$(bintype)_j"*string(j)*"p"
    fname *= ifelse(p==1,"+","-")*"t"*string(t)
    fname *= "_Nmax"*string(Nmax)*".bin"
    return fname
end

function check_binexists(JPT,params,bintype;dirname="./A3files/")
    fname = make_binname(JPT,params,bintype;dirname=dirname)
    return isfile(fname)
end

"""
cfp = <[(n12l12s12j12t12,n3j3)JT||NiJT>
"""
function get_dim_cfp!(cfpdim,cfpvals_ch,JPT,params,to)    
    j,p,t = JPT 
    nlsjts = Matrix{Int64}[ ]
    for N = 0:params.N3max
        diag = Float64[]
        for ncycle = 1:2 # 1 is to calc nphys,north, 2 is to fill nlsjtNLJ
            ndim = 0
            for n12 = 0:div(N,2)
                for n3 = 0:div(N-2*n12,2)
                    for l12 = 0:N-2*n12 -2*n3
                        l3 = N- 2*n12 -2*n3 -l12
                        if (-1)^(l12+l3) != p;continue;end
                        for s12 = 0:1
                            for j12 = abs(l12-s12):l12+s12
                                for t12 = 0:1
                                    if !tri_check(2*t12,1,t);continue;end
                                    if (-1)^(l12+s12+t12) != -1 ;continue;end
                                    for j3 = abs(2*l3-1):2:2*l3 +1
                                        if !tri_check(2*j12, j3, j);continue;end
                                        ndim += 1
                                        @assert N == 2*n12 + l12 + 2*n3 + l3 "N must be 2*n12+l12+2*n3+l3"
                                        if ncycle == 1                                            
                                            r = anti_op_isospin(params,n12,l12,s12,j12,t12,n3,l3,j3,
                                                                n12,l12,s12,j12,t12,n3,l3,j3,j,t,to)
                                            push!(diag,r)
                                        end
                                        if ncycle == 2
                                            nlsjt = nlsjts[N+1]
                                            nlsjt[1,ndim] = n12; nlsjt[2,ndim] = l12
                                            nlsjt[3,ndim] = s12; nlsjt[4,ndim] = j12
                                            nlsjt[5,ndim] = t12; nlsjt[6,ndim] = n3
                                            nlsjt[7,ndim] = l3;  nlsjt[8,ndim] = j3  
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if ncycle == 1
                push!(nlsjts,zeros(Int64,8,ndim))
                north = Int(round(sum(diag)))
                push!(cfpdim,[ndim,north])
                push!(cfpvals_ch,zeros(Float64,north,ndim))
            end
        end
    end
    return nlsjts
end

function set_cfps(cfpdim,cfpvals_ch,JPT,nlsjts,params,to)
    fname = make_binname(JPT,params,"cfp")
    io = open(fname,"w")
    j,p,t= JPT; Nmax=params.N3max
    write(io, j);write(io,p);write(io,t);write(io,Nmax)    
    j,p,t = JPT 
    for N = 0:params.N3max   
        nphys, north = cfpdim[N+1]
        nlsjt = @views nlsjts[N+1]
        if nphys == 0;continue;end
        if north == 0;continue;end
        #println("N = $N JPT $JPT i $i ndim/north $nphys $north ")
        A = zeros(Float64,nphys,nphys)   
        for ib = 1:nphys 
            n12,l12,s12,j12,t12,n3,l3,j3 = @view nlsjt[:,ib]
            for ik = 1:ib
                n45,l45,s45,j45,t45,n6,l6,j6= @view nlsjt[:,ik]
                @timeit to "anti" anti = anti_op_isospin(params,
                                                         n12,l12,s12,j12,t12,n3,l3,j3,
                                                         n45,l45,s45,j45,t45,n6,l6,j6,j,t,to)
                A[ib,ik] = A[ik,ib] = anti
            end
        end     
        @timeit to "eigen" vals,vecs = eigen(A)
        write_cfp_bin(io,north,nphys,vals,vecs,JPT,N,params,cfpvals_ch[N+1])     
    end
    close(io)
    return nothing
end

function write_cfp_bin(io,north,nphys,vals,vecs,JPT,N,params,cfpvals_ch)
    write(io, nphys);write(io,north);write(io,N)
    cfp = zeros(Float64,north,nphys)
    hit = 0
    @assert nphys == length(vals)
    for i in eachindex(vals)
        val = vals[i]
        vec = @view vecs[:,i]
        if abs(val - 1.0) < 1.e-5
            hit +=1
            cfp[hit,:] .= vec
            cfpvals_ch[hit,:] .= vec
            write(io,vec)
        elseif abs(val-0.0) < 1.e-5
            nothing
        else            
            @error "warn! something is wrong:eval cfp $val"
            exit()
        end
    end
    @assert north == hit "warn! something is wrong for north. north $north must be $hit (# of states with eval=1)"
    return nothing
end

function read_cfp_bin!(cfpdim,cfpval_ch,JPT,params,to)
    fname = make_binname(JPT,params,"cfp")
    io = open(fname,"r")
    J = read(io,Int); P = read(io,Int); T = read(io,Int); Nmax = read(io,Int)
    for N = 0:Nmax
        nphys, north = cfpdim[N+1]
        if north == 0; continue;end
        tnphys = read(io,Int)
        tnorth = read(io,Int)
        tN = read(io,Int)
        @assert nphys == tnphys "nphys $nphys != tnphys $tnphys"
        @assert north == tnorth "north $north != tnorth $tnorth"
        @assert N == tN "N $N != tN $tN"
        cfpmat = cfpval_ch[N+1]
        for i = 1:north
            vec_read = [ read(io,Float64) for j=1:nphys]
            cfpmat[i,:] .= vec_read
        end
    end
    close(io)
    return nothing
end

function JacobiHO_3NF(LECs,params3NF,JPTs,cfps,inner_integ;verbose=false,debug=true)
    N3max = params3NF.N3max
    cfpdims = cfps.cfpdims
    nallo = nallp = 0
    dict_N_idx = Dict("phys" =>[ Dict{Int64,Int64}() for idx = 1:length(JPTs)],
                      "orth" =>[ Dict{Int64,Int64}() for idx = 1:length(JPTs)])
    dict_JPT_dim = Dict{String,Dict{Int64,Int64}}("phys"=>Dict{Int64,Int64}(),"orth"=>Dict{Int64,Int64}())
    for (idx_JPT,JPT) in enumerate(JPTs)
        J,P,T = JPT        
        notmp = nptmp = 0
        for N = 0:N3max
            nphys,north = cfpdims[idx_JPT][N+1]
            dict_N_idx["orth"][idx_JPT][N] = notmp
            dict_N_idx["phys"][idx_JPT][N] = nptmp
            notmp += north
            nptmp += nphys
        end
        dict_JPT_dim["phys"][idx_JPT] = nptmp
        dict_JPT_dim["orth"][idx_JPT] = notmp
        nallp += nptmp
        nallo += notmp
        if verbose
            println("J",@sprintf("%3i",J),"  P",@sprintf("%3i",P),"  T",@sprintf("%3i",T),"   # of ME ",@sprintf("%10i",nptmp))
        end
    end
    println("# matrix element:        |NiJT> ",@sprintf("%12i",nallo))
    println("# matrix element:|(nlsjtNJT)JT> ",@sprintf("%12i",nallp))

    for (idx_JPT,JPT) in enumerate(JPTs)       
        if !check_binexists(JPT,params3NF,"JacobiHO") || debug
            calc_write_JacobiHO_ME(idx_JPT,N3max,JPT,cfps,dict_JPT_dim,dict_N_idx,params3NF,LECs,inner_integ;debug=debug)
        end
    end
    return nothing
end

function read_JacobiHO_3NF(params3NF,JPTs;verbose=false)
    ME_JacobiHO = Matrix{Float64}[ ]
    for (idx_JPT,JPT) in enumerate(JPTs)
        J,P,T = JPT
        fn = make_binname(JPT,params3NF,"JacobiHO_"*params3NF.regulator)
        io = open(fn,"r")
        dim_orth = read(io,Int64)
        Vmat = zeros(Float64,dim_orth,dim_orth)
        for i = 1:dim_orth            
            tmp = @view Vmat[i,:]
            tmp .= [read(io,Float64) for j = 1:dim_orth]
        end
        close(io)
        push!(ME_JacobiHO,Vmat)
        if verbose
            tnorm = norm(Vmat,2)
            println("   J =",@sprintf("%3i", J),"   P =",@sprintf("%3i", P),"   T =",@sprintf("%3i", T),
                    "   dim_orth ",@sprintf("%5i",dim_orth),"  norm ",@sprintf("%12.4e",tnorm))
        end
    end
    return ME_JacobiHO
end

function calc_write_JacobiHO_ME(idx_JPT::Int,N3max::Int64,JPT,
                                cfps,dict_JPT_dim,dict_N_idx,
                                params3NF,LECs,inner_integ;debug=false)
    fn = make_binname(JPT,params3NF,"JacobiHO_"*params3NF.regulator)
    J,P,T = JPT
    chs = cfps.chs
    cfpdims = cfps.cfpdims
    cfpvals = cfps.cfpvals
    dim_phys = dict_JPT_dim["phys"][idx_JPT]
    dim_orth = dict_JPT_dim["orth"][idx_JPT]
    Vmat = zeros(Float64,dim_phys,dim_phys)
    Cmat = zeros(Float64,dim_orth,dim_phys)
    tmat = zeros(Float64,dim_orth,dim_phys)
    Vret = zeros(Float64,dim_orth,dim_orth)    
    
    @threads for N_bra = 0:N3max
        chmat_bra = chs[idx_JPT][N_bra+1]
        nphys_bra,north_bra = cfpdims[idx_JPT][N_bra+1]
        cfpval_bra = cfpvals[idx_JPT][N_bra+1]
        if nphys_bra * north_bra ==0;continue;end
        ibra_p_ini = dict_N_idx["phys"][idx_JPT][N_bra] 
        ibra_p_fin = ibra_p_ini + nphys_bra
        ibra_o_ini = dict_N_idx["orth"][idx_JPT][N_bra] 
        ibra_o_fin = ibra_o_ini + north_bra
        target_C = @view Cmat[ibra_o_ini+1:ibra_o_fin,ibra_p_ini+1:ibra_p_fin]
        target_C .= cfpval_bra
        for ibra_nlsjt = 1:nphys_bra
            b_n12,b_l12,b_s12,b_j12,b_t12,b_n3,b_l3,b_j3 = @view chmat_bra[:,ibra_nlsjt]
            bra = get_ket_JacobiHO(b_n12,b_l12,b_s12,b_j12,b_t12,b_n3,b_l3,b_j3,J,T)
            idx_bra = ibra_p_ini + ibra_nlsjt
            for N_ket = N_bra:N3max
                chmat_ket = chs[idx_JPT][N_ket+1]
                nphys_ket,north_ket = cfpdims[idx_JPT][N_ket+1]
                if nphys_ket * north_ket ==0;continue;end
                iket_p_ini = dict_N_idx["phys"][idx_JPT][N_ket]
                iket_min = ifelse(N_bra==N_ket,ibra_nlsjt,1)
                for iket_nlsjt = iket_min:nphys_ket
                    k_n12,k_l12,k_s12,k_j12,k_t12,k_n3,k_l3,k_j3 = @view chmat_ket[:,iket_nlsjt] 
                    ket = get_ket_JacobiHO(k_n12,k_l12,k_s12,k_j12,k_t12,k_n3,k_l3,k_j3,J,T)
                    idx_ket = iket_p_ini + iket_nlsjt
                    v_134 = JacobiHO_3NF_TPE(params3NF,LECs,inner_integ,bra,ket)
                    v_d = JacobiHO_3NF_OPE(params3NF,LECs,inner_integ,bra,ket)
                    v_e = JacobiHO_3NF_Contact(params3NF,LECs,inner_integ,bra,ket)
                    vmat = v_134 + v_d + v_e
                  
                    Vmat[idx_bra,idx_ket] = Vmat[idx_ket,idx_bra] =  vmat
                end
            end            
        end
    end
    BLAS.gemm!('N','N',3.0,Cmat,Vmat,0.0,tmat)
    BLAS.gemm!('N','T',1.0,tmat,Cmat,0.0,Vret)

    norm_C = norm(Cmat,2)
    norm_V = norm(3.0 .* Vmat,2)
    norm_CVC = norm(Vret,2)
                
    println("   J =",@sprintf("%3i", J),"   P =",@sprintf("%3i", P),"   T =",@sprintf("%3i", T),
            "  dimp",@sprintf("%5i",dim_phys),"  dimo",@sprintf("%5i",dim_orth),
            "  normC ",@sprintf("%12.5f",norm_C),"  normV ",@sprintf("%12.5f",norm_V),
            "  normCVC ",@sprintf("%12.5f",norm_CVC))

    io = open(fn,"w")
    write(io,dim_orth)
    for i = 1:dim_orth
        write(io,Vret[:,i])
    end
    close(io)

    return nothing
end

function JacobiHO_3NF_TPE(params3NF,LECs,inner_integ,bra,ket)    
    if params3NF.regulator == "local"        
        return tpe_JacobiHO_local(params3NF,LECs,inner_integ,bra,ket)
    elseif params3NF.regulator == "nonlocal"
        @error "tpe_JacobiHO_nonlocal not defined now!!"
        #return tpe_JacobiHO_nonlocal(params3NF,LECs,inner_integ,bra,ket)
    else
        @error "regulator = $(params3NF.regulator) is not supported now"
    end    
end

function JacobiHO_3NF_OPE(params3NF,LECs,inner_integ,bra,ket)
    cD = LECs.dLECs["cD"] 
    coeffD = cD * hc * (hc/Fpi)^4 * (hc/params3NF.LambdaChi)
    
    if params3NF.regulator == "local"        
        return ope_JacobiHO_local(params3NF,coeffD,inner_integ,bra,ket)
    elseif params3NF.regulator == "nonlocal"
        return ope_JacobiHO_nonlocal(params3NF,coeffD,inner_integ,bra,ket)
    else
        @error "regulator = $(params3NF.regulator) is not supported now"
    end    
end

function JacobiHO_3NF_Contact(params3NF,LECs,inner_integ,bra,ket)
    cE = LECs.dLECs["cE"]
    coeffE = cE * hc * (hc/Fpi)^4 * (hc/params3NF.LambdaChi)
    if params3NF.regulator == "local"
        return contact_JacobiHO_local(params3NF,coeffE,inner_integ,bra,ket)
    elseif params3NF.regulator == "nonlocal"
        return contact_JacobiHO_nonlocal(params3NF,coeffE,inner_integ,bra,ket)
    else
        @error "regulator = $(params3NF.regulator) is not supported now"
    end
    return nothing
end

"""
Eq.(14) in P.Navratil
"""
function contact_JacobiHO_nonlocal(params3NF,coeffE,inner_integ::inner_integ_nonlocal,bra,ket)
    b_n12=bra.n; b_l12=bra.l; b_s12=bra.s; b_j12=bra.j; b_t12=bra.t; b_n3=bra.N; b_l3=bra.L; b_j3 = bra.J
    dJ123 = bra.dJ3; dT123 = bra.dT3
    k_n12=ket.n; k_l12=ket.l; k_s12=ket.s; k_j12=ket.j; k_t12=ket.t; k_n3=ket.N; k_l3=ket.L; k_j3 = ket.J
    if dJ123 != 1; return 0.0;end
    if k_l12 != 0; return 0.0;end
    if b_l12 != 0; return 0.0;end
    if k_l3 != 0; return 0.0;end
    if b_l3 != 0; return 0.0;end
    if b_s12 != b_j12; return 0.0;end
    if k_s12 != k_j12; return 0.0;end
    if b_s12 != k_s12; return 0.0;end
    if b_t12 == k_t12 == 0; return 0.0; end
    tdot = IsospinDot(b_t12,k_t12,dT123,params3NF.dWS.d6j_lj)
    if tdot==0.0; return 0.0; end
    e_bra = 2*b_n12 + 2*b_n3; tkey1 = get_nkey3(e_bra,b_n12,b_n3); freg1 = inner_integ.integ_E[tkey1]
    e_ket = 2*k_n12 + 2*k_n3; tkey2 = get_nkey3(e_ket,k_n12,k_n3); freg2 = inner_integ.integ_E[tkey2]
    v_e = coeffE / (12.0 *sqrt(3.0) *pi^4)* tdot * freg1 * freg2
    return v_e
end

function prep_Z0s_Z0Xs_fKs(params)
    N3max = params.N3max
    Xmax = 2*N3max # will be truncated
    ps = params.meshpoints.ps
    ws = params.meshpoints.wps
    Z0s  = Dict{Int64,Float64}() 
    Z0Xs = Dict{Int64,Float64}()
    fKs  = Dict{Int64,Float64}() 
    fKXs  = Dict{Int64,Float64}() 

    for (imesh1,w1) in enumerate(ws)
        p = ps[imesh1]
        xi1 = sqrt(2.0) * p
        Z0s[imesh1]= Z0_local_reg(xi1,params.Lambda3NF,ps,ws)
        for K = 0:1:2
            fKs[get_nkey2(K,imesh1)] = f_K(K,xi1,params.Lambda3NF,ps,ws)
        end
        xi1 = p /sqrt(2.0)
        for (imesh2,w2) in enumerate(ws)
            xi2 = ps[imesh2] * sqrt(3/2)
            for X = 0:Xmax
                meshkey = get_nkey3(X,imesh1,imesh2)               
                Z0Xs[meshkey]=Z0X_local_reg(X,xi1,xi2,params.Lambda3NF,ps,ws)
            end
            for K = 0:1:2
                for X = 0:Xmax+1
                    tmp = Float64(f_KX(K,X,xi1,xi2,params))
                    fKXs[get_nkey4(K,X,imesh1,imesh2)] = tmp  
                end
            end
        end

    end
    return Z0s,Z0Xs,fKs,fKXs
end

"""
``Z_0(r;\\Lambda)=\\frac{1}{2\\pi^2}\\int dq q^2 j_0(qr)F(q^2;\\Lambda)`` Eq.(13) in Few Body Syst. (2007) 41:117-140
"""
function Z0_local_reg(r,Lambda,qs,ws;npower=2,debug=false)
    nmesh = length(qs)
    qinteg = 0.0
    for i = 1:nmesh
        q = qs[i]
        w = ws[i]
        F = exp( - (q^(2*npower) / (Lambda/hc)^4 ) )
        qinteg += w * q^2 * sphericalbesselj(0,q*r) * F
        if debug
            println("i",@sprintf("%5i",i),"   w",@sprintf("%12.5f",w),
            "   q2 ",@sprintf("%12.5f",q^2),"   F ",@sprintf("%12.5f",F), 
            "   j0 ",@sprintf("%12.5f",sphericalbesselj(0,q*r)))
        end
    end
    Z = 1.0/(2*pi^2) * qinteg
    return Z 
end

"""
``Z_{0,X}(r_1,r_2;\\Lambda) = \\frac{1}{2\\pi^2} \\int dq q^2 j_X(qr_1) j_X(qr_2) F(q^2;\\Lambda) ``
Eq.(16) in Few Body Syst. (2007) 41:117-140.
"""
function Z0X_local_reg(X,r1,r2,Lambda,qs,ws;npower=2,debug=false)
    nmesh = length(qs)
    qinteg = 0.0
    for i = 1:nmesh
        q = qs[i]
        w = ws[i]
        F = exp( - (q^(2*npower)/(Lambda/hc)^4) )
        qinteg += w * q^2 * sphericalbesselj(X,q*r1) * sphericalbesselj(X,q*r2) * F
        sqj = sphericalbesselj(X,q*r1) * sphericalbesselj(X,q*r2) 
        if debug
            println("X $X r1 $r1 r2 $r2 jXsq $sqj F $F")
        end
    end
    Z0X = 1.0/(2*pi^2) * qinteg
    if debug
        println("Z0X $Z0X")
    end
    return Z0X
end

"""
Navratil Eq.(39):  

``
f_K(r ; \\Lambda)=\\frac{1}{2 \\pi^2} \\int d q q^2 j_K(q r) \\frac{q^2 F\\left(q^2 ; \\Lambda\\right)}{q^2+M_\\pi^2}
``
"""
function f_K(K,r,Lambda,qs,ws;npower=2,debug=false)
    Mpi_fm = mean(mpis) / hc
    qinteg = 0.0
    for (i,w) in enumerate(ws)
        q = qs[i]
        q2= q^2
        F = exp( - (q2^(npower) / (Lambda/hc)^4 ) )
        qinteg += w * q2 * sphericalbesselj(K,q*r) * F * ifelse(K==1,q,q2) / (q2+Mpi_fm^2)
        if debug
            println("i",@sprintf("%5i",i)," K",@sprintf("%5i",K),"   w",@sprintf("%12.5f",w),
            "   q2 ",@sprintf("%12.5f",q^2),"   F ",@sprintf("%12.5f",F), 
            "   jK ",@sprintf("%12.5f",sphericalbesselj(K,q*r)))
        end
    end
    Z = 1.0/(2*pi^2) * qinteg
    return Z 
end

function f_KX(K,X,xi1,xi2,params;debug=false)
    Lambda = params.Lambda3NF
    qs = params.meshpoints.ps
    ws = params.meshpoints.wps
    us = params.meshpoints.us
    wus = params.meshpoints.wus
    uinteg = 0.0
    for (i,wu) in enumerate(wus)
        u = us[i]    
        r_in = sqrt(xi1^2+xi2^2-2*xi1*xi2*u)
        nume = f_K(K,r_in,Lambda,qs,ws)        
        deno = r_in^K
        uinteg += 0.5 * wu * Legendre(X,u) *nume/deno
        if debug
        end
    end
    return uinteg
end

"""
Eq.(15) in P.Navratil
"""
 function contact_JacobiHO_local(params3NF,coeffE,inner_integ::inner_integ_local,bra::ket_JacobiHO,ket::ket_JacobiHO)
    b_n12=bra.n; b_l12=bra.l; b_s12=bra.s; b_j12=bra.j; b_t12=bra.t; b_n3=bra.N; b_l3=bra.L; b_j3 = bra.J
    dJ123 = bra.dJ3; dT123 = bra.dT3
    k_n12=ket.n; k_l12=ket.l; k_s12=ket.s; k_j12=ket.j; k_t12=ket.t; k_n3=ket.N; k_l3=ket.L; k_j3 = ket.J
    if b_s12 != k_s12; return 0.0; end
    if b_t12 == k_t12 == 0; return 0.0; end
    tdot = IsospinDot(b_t12,k_t12,dT123,params3NF.dWS.d6j_lj)
    d6j_int = params3NF.dWS.d6j_int
    d6j_lj = params3NF.dWS.d6j_lj
    dcgm0 = params3NF.dWS.dcgm0
    hats  = hat(b_j12)*hat(k_j12)*hat(b_j3/2)*hat(k_j3/2)
    phase = (-1)^(div(dT123-1,2)+div(abs(b_j3-k_j3),2)+b_l12+b_l3+b_s12)
    
    Xmin = max(max(abs(b_l12-k_l12),abs(b_j12-k_j12)),div(abs(b_j3-k_j3),2))
    Xmax = min(min(b_l12+k_l12,b_j12+k_j12),b_j3+k_j3)

    e1 = bra.e1 ; e2 = bra.e2; ebra = bra.E
    e1p= ket.e1 ; e2p= ket.e2; eket = ket.E

    Ekey = get_nkey2(ebra,eket)
    e4key = get_nkey4(e1,e2,e1p,e2p)
    key_bra = get_nkey4(b_n12,b_l12,b_n3,b_l3)
    key_ket = get_nkey4(k_n12,k_l12,k_n3,k_l3)
    integ_dict = inner_integ.integ_E[Ekey][e4key][key_bra][key_ket]

    X6js = 0.0
    for X = Xmin:Xmax
        if !tri_check(b_l12,k_l12,X);continue;end
        if !tri_check(b_j12,k_j12,X);continue;end
        if !tri_check(b_j3,k_j3,2*X);continue;end
        if !tri_check(b_l3,k_l3,X);continue;end

        x6j_1 = call_d6j_nond(k_l12,b_l12,X,b_j12,k_j12,b_s12,d6j_int)
        x6j_2 = call_d6j(b_j12*2,k_j12*2,X*2,k_j3,b_j3,dJ123,d6j_lj)
        x6j_3 = call_d6j(b_l3*2,k_l3*2,X*2,k_j3,b_j3,1,d6j_lj)
        cgs = call_dcgm0(k_l12,X,b_l12,dcgm0) *call_dcgm0(k_l3,X,b_l3,dcgm0) * hat(k_l12)*hat(k_l3) 
        
        tsum = (-1)^X * hat(X)^2 * x6j_1 * x6j_2* x6j_3 * cgs * integ_dict[X]
        X6js += tsum
    end
    return coeffE * tdot * (hats * phase * X6js)
end

function flip_needed_6j(j1,j2,J12,j3,J,J23,dict)
    ret = 0.0
    if j1 <= j2
        key2_1 = get_nkey3(j1,j2,J12)
        key2_2 = get_nkey3(j3,J,J23)
        ret = dict[key2_1][key2_2]
    else
        key2_1 = get_nkey3(j2,j1,J12)
        key2_2 = get_nkey3(J,j3,J23)
        ret = dict[key2_1][key2_2]
    end
    return ret
end

"""
Navratil Eq.(30)
(-1)^(-t12-t12') is needed to cancel the factor from IsospinDot
"""
function ope_JacobiHO_nonlocal(params3NF,coeffD,inner_integ::inner_integ_nonlocal,bra::ket_JacobiHO,ket::ket_JacobiHO)
    b_n12=bra.n; b_l12=bra.l; b_j12=bra.j; b_t12=bra.t; b_n3=bra.N; b_l3=bra.L; b_j3 = bra.J
    dJ123 = bra.dJ3; dT123 = bra.dT3
    k_n12=ket.n; k_l12=ket.l; k_j12=ket.j; k_t12=ket.t; k_n3=ket.N; k_l3=ket.L; k_j3 = ket.J
    if k_l12 != 0; return 0.0;end
    if b_l12 != 0; return 0.0;end
    if abs(b_l3-k_l3) > 2; return 0.0; end
    if !tri_check(b_j12,k_j12,1); return 0.0;end
    tdot = IsospinDot(b_t12,k_t12,dT123,params3NF.dWS.d6j_lj)
    if tdot == 0.0; return 0.0; end

    d6j_int = params3NF.dWS.d6j_int
    d6j_lj = params3NF.dWS.d6j_lj
    d9j_lsj = params3NF.dWS.d9j_lsj
    dcgm0 = params3NF.dWS.dcgm0

    hats = hat(b_j12)*hat(k_j12)*hat(b_j3/2)*hat(k_j3/2)
    phase = (-1)^(b_n12+b_n3+k_n12+k_n3 + (dJ123-b_j3+b_l3+k_l3)/2 + b_l3 - b_t12 - k_t12)

    w6js  = call_d6j(b_j12*2,k_j12*2,1*2,1,1,1,d6j_lj) 
    w6js *= call_d6j(b_j12*2,k_j12*2,1*2,k_j3,b_j3,dJ123,d6j_lj)
    if w6js == 0.0; println("w6js0");return 0.0;end
    
    e1 = bra.e1 ; e2 = bra.e2; Ebra = bra.E
    e1p= ket.e1 ; e2p= ket.e2; Eket = ket.E
  
    Ekey = get_nkey2(Ebra,Eket)
    e4key = get_nkey4(e1,e2,e1p,e2p)
    keybra = get_nkey4(b_n12,b_l12,b_n3,b_l3)
    keyket = get_nkey4(k_n12,k_l12,k_n3,k_l3)
    dict_integ = inner_integ.integ_D[Ekey][e4key][keybra][keyket]

    ret = 0.0
    for K = 0:2:2
        if !tri_check(b_l3,k_l3,K);continue;end
        K9j = call_d9j_lsj(K*2,1*2,1*2,b_l3*2,1,b_j3,k_l3*2,1,k_j3,d9j_lsj)
        Kfac = hat(K) * call_dcgm0(1,1,K,dcgm0) * K9j
        for K1=0:K            
            K1fac =hat(K-K1) * sqrt(binomial(2*K+1,2*K1))*  (-1)^(b_l3+K1) 
            Xmin = max(abs(b_l3-K1),abs(k_l3-(K-K1)))
            Xmax = min(b_l3+K1,k_l3+K-K1)
            for X = Xmin:Xmax
                Xfac  = hat(X)* hat(k_l3) * call_dcgm0(K1,X,b_l3,dcgm0) * call_dcgm0(k_l3,K-K1,X,dcgm0) # *dcgm0[get_nkey3(K1,X,b_l3)] * dcgm0[get_nkey3(k_l3,K-K1,X)] 
                Xfac *= call_d6j_nond(k_l3,K-K1,X,K1,b_l3,K,d6j_int) #d6j_int[X+1][key1]
                integ = dict_integ[get_nkey3(K,K1,X)]
                ret += Kfac * K1fac * Xfac * integ
            end
        end
    end

    v_d = - abs(gA) * coeffD / (12.0 *sqrt(3.0) *pi^4) * tdot * hats * phase *w6js * ret
    return v_d
end

"""
Navratil Eq.(38)
"""
function ope_JacobiHO_local(params3NF,coeffD,inner_integ::inner_integ_local,bra::ket_JacobiHO,ket::ket_JacobiHO)
    b_n12=bra.n; b_l12=bra.l; b_s12=bra.s; b_j12=bra.j; b_t12=bra.t; b_n3=bra.N; b_l3=bra.L; b_j3 = bra.J
    dJ123 = bra.dJ3; dT123 = bra.dT3
    k_n12=ket.n; k_l12=ket.l; k_s12=ket.s; k_j12=ket.j; k_t12=ket.t; k_n3=ket.N; k_l3=ket.L; k_j3 = ket.J
    if b_s12 == k_s12 ==0; return 0.0;end
    tdot = IsospinDot(b_t12,k_t12,dT123,params3NF.dWS.d6j_lj)
    if tdot == 0.0; return 0.0; end

    d6j_int = params3NF.dWS.d6j_int
    d6j_lj = params3NF.dWS.d6j_lj
    d9j_int = params3NF.dWS.d9j_int
    d9j_lsj = params3NF.dWS.d9j_lsj
    dcgm0 = params3NF.dWS.dcgm0

    hats = hat(b_j12)*hat(k_j12)*hat(b_j3/2)*hat(k_j3/2)*hat(b_s12)*hat(k_s12) * hat(k_l12)*hat(k_l3)
    phase6j = (-1)^((dJ123-b_j3)/2 +b_s12+k_j12) * call_d6j(b_s12*2,k_s12*2,2,1,1,1,d6j_lj) #d6j_dot[1][get_nkey2(b_s12,k_s12)]

    e1 = bra.e1 ; e2 = bra.e2; Ebra = bra.E
    e1p= ket.e1 ; e2p= ket.e2; Eket = ket.E

    Ekey = get_nkey2(Ebra,Eket)
    e4key = get_nkey4(e1,e2,e1p,e2p)
    keybra = get_nkey4(b_n12,b_l12,b_n3,b_l3)
    keyket = get_nkey4(k_n12,k_l12,k_n3,k_l3)
    dict_integ = inner_integ.integ_D[Ekey][e4key][keybra][keyket]

    ret = 0.0
    for K = 0:2:2
        Kfac = hat(K) * (-1)^(K/2) * call_dcgm0(1,1,K,dcgm0) 
        for V = abs(b_l12-k_l12):b_l12+k_l12
            Vfac = (-1)^(V+k_l12) * hat(b_l12) *  call_dcgm0(b_l12,k_l12,V,dcgm0) 

            Xmin = max(abs(K-V),abs(b_l3-k_l3))
            Xmax = min(K+V,b_l3+k_l3)
            for X = Xmin:Xmax
                Xfac = hat(X)^2 * call_dcgm0(X,K,V,dcgm0) * call_dcgm0(X,k_l3,b_l3,dcgm0)
                Zmin = max(max(max(abs(X-1),abs(V-1)),abs(b_j12-k_j12)),abs(div(b_j3-k_j3,2)))
                Zmax = min(min(min(V+1,X+1),b_j12+k_j12),div(b_j3+k_j3,2))
                integ = dict_integ[get_nkey2(K,X)]                

                for Z = Zmin:Zmax
                    z6j1 = call_d6j(b_j12*2,k_j12*2,Z*2,k_j3,b_j3,dJ123,d6j_lj)
                    if z6j1 == 0.0;continue;end

                    z6j2 = call_d6j_nond(V,1,Z,1,X,K,d6j_int)
                    if z6j2 ==0.0; continue;end

                    z9j2 = call_d9j_lsj(X*2,1*2,Z*2,b_l3*2,1,b_j3,k_l3*2,1,k_j3,d9j_lsj)
                    z9j1 = call_d9j_int_notdoubled(V,1,Z,b_l12,b_s12,b_j12,k_l12,k_s12,k_j12,d9j_int)

                    Zfac = hat(Z)^2 * z6j1 * z6j2 * z9j1 *z9j2
                    
                    ret += Kfac * Vfac * Xfac * Zfac * integ
                end
            end
        end
    end
    v_d = -abs(gA) * coeffD * 3/2 *tdot *hats * phase6j *ret
    return v_d
end

function tpe_JacobiHO_local(params3NF,LECs,inner_integ::inner_integ_local,bra::ket_JacobiHO,ket::ket_JacobiHO;debug=false)
    b_n12=bra.n; b_l12=bra.l; b_s12=bra.s; b_j12=bra.j; b_t12=bra.t; b_n3=bra.N; b_l3=bra.L; b_j3 = bra.J
    dJ123 = bra.dJ3; dT123 = bra.dT3
    k_n12=ket.n; k_l12=ket.l; k_s12=ket.s; k_j12=ket.j; k_t12=ket.t; k_n3=ket.N; k_l3=ket.L; k_j3 = ket.J
    if b_s12 == k_s12 ==0; return 0.0;end
    tdot = IsospinDot(b_t12,k_t12,dT123,params3NF.dWS.d6j_lj)
    if tdot == 0.0; return 0.0; end
    Mpi_fm = mean(mpis) / hc

    coeff_c1 = - 6 * Mpi_fm^2 * gA2 / (Fpi/hc)^4 * LECs.dLECs["c1_NNLO"] * (hc2 /1000.0) 
    coeff_c3 = 3* gA2 / (Fpi/hc)^4 * LECs.dLECs["c3_NNLO"] * hc2 / 1000.0
    coeff_c4 = - 216  * gA2 / (4*(Fpi/hc)^4) *  LECs.dLECs["c4_NNLO"] * hc2 / 1000.0

    d6j_int = params3NF.dWS.d6j_int
    d6j_lj  = params3NF.dWS.d6j_lj

    d9j_lsj = params3NF.dWS.d9j_lsj
    d9j_int = params3NF.dWS.d9j_int
    dcgm0 = params3NF.dWS.dcgm0

    hats = hat(b_j12)*hat(k_j12)*hat(b_j3/2)*hat(k_j3/2)*hat(b_s12)*hat(k_s12) * hat(k_l12)*hat(k_l3)
    s6j = call_d6j(b_s12*2,k_s12*2,1*2,1,1,1,d6j_lj) 
    phase6j = (-1)^((dJ123-b_j3)/2 +b_s12+k_j12) 

    rac_9j = call_d9j_lsj(2,2,2,k_t12*2,1,1,b_t12*2,1,1,d9j_lsj)
    rac_9jc4= rac_9j / call_d6j(b_t12*2,k_t12*2,2,1,1,1,d6j_lj)

    e1 = bra.e1 ; e2 = bra.e2; Ebra = bra.E
    e1p= ket.e1 ; e2p= ket.e2; Eket = ket.E
    Ekey = get_nkey2(Ebra,Eket)
    e4key = get_nkey4(e1,e2,e1p,e2p)
    keybra = get_nkey4(b_n12,b_l12,b_n3,b_l3)
    keyket = get_nkey4(k_n12,k_l12,k_n3,k_l3)
    dict_c1 = inner_integ.integ_c1[Ekey][e4key][keybra][keyket]
    dict_c34= inner_integ.integ_c34[Ekey][e4key][keybra][keyket]
    ret_c1 = ret_c3 = ret_c4 = 0.0
    for V = abs(b_l12-k_l12):b_l12+k_l12
        Vfac = hat(V) * call_dcgm0(V,k_l12,b_l12,dcgm0) 
        if Vfac == 0.0; continue;end
        for R = abs(b_l3-k_l3):b_l3+k_l3
            Rfac = hat(R) * call_dcgm0(R,k_l3,b_l3,dcgm0)
            if Rfac == 0.0; continue;end
            ## c1 term
            if s6j != 0.0
                Ymin = max(max(abs(V-1),abs(b_j12-k_j12)),div(abs(b_j3-k_j3),2))
                Ymax = min(V+1,min(b_j12+k_j12,div(b_j3+k_j3,2)))
                for Y = Ymin:Ymax
                    if !tri_check(R,1,Y);continue;end
                    if !tri_check(b_j3/2,k_j3/2,Y);continue;end
                    V9j = call_d9j_int_notdoubled(V,1,Y,b_l12,b_s12,b_j12,k_l12,k_s12,k_j12,d9j_int)
                    R9j = call_d9j_lsj(R*2,1*2,Y*2,b_l3*2,1,b_j3,k_l3*2,1,k_j3,d9j_lsj)
                    y6j = call_d6j(b_j12*2,k_j12*2,Y*2,k_j3,b_j3,dJ123,d6j_lj)
                    if y6j == 0.0;continue;end
                    Yfac = (-1)^Y * hat(Y) * call_dcgm0(Y,1,V,dcgm0) * V9j * R9j *y6j 
                    if Yfac == 0.0; continue;end

                    for K3 = 0:1
                        K3fac = sqrt(binomial(3,2*K3)) * hat(1-K3)    
                        Xmin = max(abs(K3-Y),abs(R+K3-1))
                        Xmax = min(K3+Y,R+1-K3)
                        for X = Xmin:Xmax
                            Xfac  = hat(X)^2 * call_dcgm0(X,K3,Y,dcgm0) * call_dcgm0(X,1-K3,R,dcgm0)
                            Xfac *= call_d6j_nond(Y,X,K3,1-K3,1,R,d6j_int)
                            integ = dict_c1[get_nkey2(K3,X)]                
                            tsum = Vfac * Rfac * Yfac * K3fac * Xfac * integ
                            ret_c1 += tsum
                        end
                    end
                end
            end
            ## c3&c4 term
            c3fac = s6j
            c4fac = (-1)^(-b_t12 + b_j3 - b_s12) * rac_9jc4
            for K1 = 0:2:2
                for K2 = 0:2:2
                    # c3 term
                    Ymin = max(abs(K1-V),abs(K2-R))
                    Ymax = min(K1+V,K2+R)
                    for Y = Ymin:Ymax
                        if !tri_check(K1,V,Y); continue;end
                        if !tri_check(K2,R,Y); continue;end
                        K12fac = (-1)^((K1+K2)/2) * hat(K1) * hat(K2) * call_dcgm0(1,1,K1,dcgm0) * call_dcgm0(1,1,K2,dcgm0) 
                        Yfac = hat(Y) * call_dcgm0(Y,K1,V,dcgm0) 
                        Zmin = max(abs(Y-1),abs(V-1),abs(R-1),abs(b_j12-k_j12),div(abs(b_j3-k_j3),2))
                        Zmax = min(    Y+1,     V+1,     R+1,     b_j12+k_j12,div(b_j3+k_j3,2))
                        for Z = Zmin:Zmax
                            Zfac = hat(Z)^2 * (-1)^Z 
                            R9j = call_d9j_lsj(R*2,1*2,Z*2,b_l3*2,1,b_j3,k_l3*2,1,k_j3,d9j_lsj)
        
                            Z6j_1 = call_d6j(b_j12*2,k_j12*2,Z*2,k_j3,b_j3,dJ123,d6j_lj)
                            if Z6j_1 == 0.0;continue;end
                            Z6j_2 = call_d6j_nond(Z,1,Y,K2,R,1,d6j_int)
                            zc3_6j = call_d6j_nond(Z,1,Y,K1,V,1,d6j_int)
                            zc3_9j = call_d9j_int_notdoubled(V,1,Z,b_l12,b_s12,b_j12,k_l12,k_s12,k_j12,d9j_int)
                             
                            Z_c4 = Zfac * R9j * Z6j_1 * Z6j_2
                            Z_c3 = Z_c4 * zc3_9j * zc3_6j

                            for K3 = 0:K2
                                K3fac = sqrt(binomial(2*K2+1,2*K3)) * hat(K2-K3) 

                                Xmin = max(abs(K3-Y),abs(R-K2+K3))
                                Xmax = min(K3+Y,R+K2-K3)
                                for X = Xmin:Xmax
                                    if !tri_check(X,K3,Y);continue;end
                                    if !tri_check(X,K2-K3,R);continue;end
                                    Xfac  = hat(X)^2 * call_dcgm0(X,K3,Y,dcgm0) * call_dcgm0(X,K2-K3,R,dcgm0) 
                                    Xfac *= call_d6j_nond(Y,K2,R,K2-K3,X,K3,d6j_int)    
                                    if Xfac == 0.0; continue;end                               
                                    integ = dict_c34[get_nkey4(K1,K2,K3,X)]
                                    ret_c3 += c3fac * K12fac * Vfac * Rfac * Yfac * Z_c3 * K3fac * Xfac * integ

                                    for K4 = 0:2
                                        if !tri_check(K4,b_s12,k_s12);continue;end
                                        if !tri_check(K4,V,Z);continue;end
                                        if !tri_check(K4,K1,1);continue;end
                                        K4fac = hat(K4)^2 
                                        K4fac *= call_d9j_int_notdoubled(V,K4,Z,b_l12,b_s12,b_j12,k_l12,k_s12,k_j12,d9j_int)
                                        K4fac *= call_d9j_lsj(K4*2,1*2,1*2,b_s12*2,1,1,k_s12*2,1,1,d9j_lsj)
                                        K4fac *= call_d6j_nond(Z,1,Y,K1,V,K4,d6j_int)
                                        K4fac *= call_d6j_nond(1,1,1,1,K4,K1,d6j_int)
                                        ret_c4 += c4fac * K12fac * Vfac * Rfac * Yfac * Xfac* Z_c4 * K3fac * K4fac* integ
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    v_c1 = coeff_c1 * tdot *hats * phase6j * s6j * ret_c1
    v_c3 = coeff_c3 * tdot *hats * phase6j * ret_c3
    v_c4 = coeff_c4 * tdot *hats * phase6j * ret_c4
    return v_c1 + v_c3 + v_c4
end
