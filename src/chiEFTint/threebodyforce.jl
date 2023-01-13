"""
    struct `nine_j_stuff` 9j dict and keys used for HOBs in 3NF
"""
struct nine_j_stuff
    dict9j_HOB::Vector{Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}}
    tkeys9j::Vector{Vector{Int64}}
end

struct dWS3n
    dcgm0::Dict{Int64,Float64}
    d6j_hfspin::Dict{Int64,Float64}
    d6j_int::Vector{Dict{Int64,Float64}}
    d6j_dot::Dict{Int64,Dict{Int64,Float64}}
    d6j_lj3::Dict{Int64,Dict{Int64,Float64}}
    d9j::Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}
    d9j_int::Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}
    dtri::Dict{Vector{Int64},Float64}
    keycg::Vector{Vector{Int64}}
    key6j::Vector{Int64}
    key9j::Vector{Int64}
    ninej_stuff::nine_j_stuff
end 

struct mesh_3NF
    ps::Vector{Float64}
    wps::Vector{Float64}
    qs::Vector{Float64}
    wqs::Vector{Float64}
    pqs::Matrix{Float64}
    wpqs::Vector{Float64}
    us::Vector{Float64}
    wus::Vector{Float64}
end

struct params3b    
    e3max::Int64
    N3max::Int64
    L3max::Int64
    j3max::Int64
    J12max::Int64
    hw::Int64
    chiorder3NF::Int64
    Lambda3NF::Float64
    LambdaChi::Float64
    TFrenorm::Bool
    pmax3::Float64
    rmax3::Float64
    meshpoints::mesh_3NF
    LECs::LECs
    dWS::dWS3n
    regulator::String
end

struct channel3b_Jj
    dJ123::Int64
    dT123::Int64
    P123::Int64
    L12::Int64
    S12::Int64
    J12::Int64
    T12::Int64
    dj3::Int64
    l3::Int64
end

struct state_Jacobi
    p::Float64
    q::Float64
    alpha::channel3b_Jj
end

struct ket_JacobiHO
    n::Int64
    l::Int64
    s::Int64
    j::Int64
    t::Int64
    N::Int64
    L::Int64
    J::Int64
    dJ3::Int64
    dT3::Int64
    e1::Int64
    e2::Int64
    E::Int64
end

function get_ket_JacobiHO(n,l,s,j,t,N,L,J,dJ3,dT3)
    e1 = 2*n+l
    e2 = 2*N+L
    Eket = e1+e2
    return ket_JacobiHO(n,l,s,j,t,N,L,J,dJ3,dT3,e1,e2,Eket)
end

function test3NF(;target_LECs=["c1_NNLO","c3_NNLO","c4_NNLO","cD","cE"],fn_params="optional_parameters.jl")
    to = TimerOutput()
    io =  select_io(false,"",[])
    paramsNN = init_chiEFTparams(;io=io,fn_params=fn_params)
    LECs = read_LECs(paramsNN.pottype)

    ### parameters  (to be included in params)
    TFrenorm = false

    main_basis="Jacobi" #K.Hebeler    
    main_basis="JacobiHO" #P.Navratil
  
    regulator = "nonlocal"
    regulator = "local"

    Np = Nq = 10 #25
    Nangle = 10 #25
    pmax3 = 10.0
    rmax3 = 10.0
    e3max = 8
    L3max = 0 # (notused now)
    N3max = 8 # to be 40?
    j3max = 9 # said to be sufficient, but can be e3max*2 + 3
    J12max = 7    
    chiorder3b = 2 # NNLO
    Lambda3NF = 400.0
    LambdaChi = 700.0

    @timeit to "precalc WignerSymbols" dWS = prep_dWS3n(paramsNN,e3max,N3max,J12max,j3max,to)
    meshpoints = make_3b_mesh_mats(pmax3,Np,Nq,Nangle)
    params3N = params3b(e3max,N3max,L3max,j3max,J12max,hw,chiorder3b,Lambda3NF,LambdaChi,TFrenorm,pmax3,rmax3,meshpoints,LECs,dWS,regulator)

    if main_basis == "Jacobi" # K.Hebeler type                
        Calc_3NF_in_Jacobi_coordinate(params3N,LECs,to)
    elseif main_basis == "JacobiHO" # P.Navratil type
        Calc_3NF_in_JacobiHO_coorrdinate(params3N,LECs,to)
    else
        @error "main_basis=$main_basis is not supported"
    end

    show(to, allocations = true,compact = false);println("")
end

function make_3b_mesh_mats(pmax3,Np,Nq,Nangle)
    ps,wps = Gauss_Legendre(0.0,pmax3,Np)
    qs,wqs = Gauss_Legendre(0.0,pmax3,Nq)
    pqs,wpqs,dim = make_pq(ps,wps,qs,wqs)
    us,wus = Gauss_Legendre(-1.0,1.0,Nangle)
    return mesh_3NF(ps,wps,qs,wqs,pqs,wpqs,us,wus)
end

"""
    prep_dWS3n(e3max,N3max)

    to prepare Wigner symbols (CG/6j/9j) for three body force
    Note:
    * only j1 <= j2 case is considered for 6j symbols {j1 j2 j12;j3 J j23}
    * PrecalcHOB is called in this function
"""
function prep_dWS3n(paramsNN,e3max,N3max,J12max,j3max,to)    
    dcgm0 = Dict{Int64,Float64}()
    dtri = Dict{Vector{Int64},Float64}()
    d6j_hfspin = Dict{Int64,Float64}()
    d9j = Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}()
    d9j_int = Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}()
    @timeit to "trinomial&cgm0" begin
        #(l1,0,l2,0|l,0) type CG coeff., which appear in HOBs and PW 3NF matrix elements
        lmax = 2*N3max
        for l = 0:lmax
            for l1 = 0:lmax
                for l2 = 0:lmax
                    key = [l1,l2,l]
                    dtri[key] = trinomial(l1,l2,l)
                end
            end
        end
        for l = 0:lmax
            for l1 = 0:lmax
                for l2 = abs(l-l1):l+l1
                    if !tri_check(l1,l2,l);continue;end
                    tkey = get_nkey3(l1,l2,l)
                    dcgm0[tkey] = clebschgordan(Float64,l1,0,l2,0,l,0)
                end
            end
        end
    end
    # coupling of spin 1/2 three particle (used for S12,S123,T12,T123)
    s3 = 1//2
    for s12=0:1
        for s45=0:1
            smin = max(abs(s12-s3),abs(s45-s3))
            smax = min(s12+s3,s45+s3)
            for S = smin:smax
                key = get_nkey3(s12,2*S,s45)
                d6j_hfspin[key] = wigner6j(Float64,1//2,1//2,s12,1//2,S,s45)
            end
        end
    end
    L12max = N3max*2
    @timeit to "9j" for l12 = 0:L12max
        for s12 = 0:1
            for j12 = abs(l12-s12):l12+s12
                nkey1 = Int(get_nkey3(l12,s12,j12))
                if get(d9j,nkey1,nothing) == nothing; d9j[nkey1]=Dict{Int64,Dict{Int64,Float64}}();end
                for l3 = 0:N3max
                    ds3 = 1
                    for dj3 = abs(2*l3-ds3):2:2*l3+ds3
                        j3 = dj3//2
                        nkey2 = Int(get_nkey2(l3,dj3))
                        if get(d9j[nkey1],nkey2,nothing) == nothing; d9j[nkey1][nkey2]=Dict{Int64,Float64}();end
                        for Lam = abs(l12-l3):l12+l3
                            if !tri_check(Lam,l12,l3);continue;end
                            for dS = abs(2*s12-ds3):2:2*s12+ds3
                                S = dS//2
                                if !tri_check(S,s12,s3);continue;end
                                for dJ = abs(2*j12-dj3):2:2*j12+dj3
                                    J = dJ //2
                                    if !tri_check(J,j12,j3);continue;end                                    
                                    t9j = wigner9j(l12,s12,j12,l3,s3,j3,Lam,S,J)  
                                    nkey3 = Int(get_nkey3(Lam,dS,dJ))
                                    if get(d9j[nkey1][nkey2],nkey3,nothing) == nothing; d9j[nkey1][nkey2][nkey3]=0.0;end
                                    d9j[nkey1][nkey2][nkey3] = t9j
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    @timeit to "9jint" begin 
        for l = 0:N3max
            for s = 0:1
                for j = abs(l-s):l+s
                    key1 = get_nkey3(l,s,j)
                    d9j_int[key1] = Dict{Int64,Dict{Int64,Float64}}()
                    for lp = 0:N3max
                        for sp = 0:1
                            for jp = abs(lp-sp):lp+sp
                                key2 = get_nkey3(lp,sp,jp)
                                d9j_int[key1][key2] = Dict{Int64,Float64}()                                
                                for V = abs(l-lp):l+lp
                                    for Z = abs(V-1):V+1
                                        for S = abs(s-sp):s+sp
                                            key3 = get_nkey3(V,S,Z)
                                            t9j = wigner9j(l,s,j,lp,sp,jp,V,S,Z)
                                            d9j_int[key1][key2][key3] = t9j
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

    @timeit to "6jsymbols" begin
        Jmax = 2*N3max
        int_jmax = 2*N3max 
        d6j_int = [ Dict{Int64,Float64}() for J = 0:Jmax]
        d6j_lj3  = Dict{Int64,Dict{Int64,Float64}}()
        for totJ = 0:Jmax
            tdict = d6j_int[totJ+1]
            for j1 = 0:int_jmax
                for j2 = j1:int_jmax
                    if !tri_check(j1,j2,totJ);continue;end
                    # integer-only part
                    for Jp = 0:Jmax
                        for j3 = 0:int_jmax
                            if !tri_check(j2,j3,Jp);continue;end
                            for j4 = 0:int_jmax
                                if !tri_check(j3,j4,totJ);continue;end                        
                                if !tri_check(j1,j4,Jp);continue;end
                                t6j = wigner6j(Float64,j1,j2,totJ,j3,j4,Jp)
                                nkey = get_nkey_from_key6j(j1,j2,j3,j4,Jp)
                                tdict[nkey] = t6j
                            end
                        end
                    end
                    # 3-integer + 3-half-integer part{j1,j2,totJ; J3',J3,J123}
                    key_int = get_nkey3(j1,j2,totJ)
                    d6j_lj3[key_int] = Dict{Int64,Float64}()
                    for J3p = 1:2:Jmax*2
                        for J3 = 1:2:Jmax*2
                            if !tri_check(J3p//2,J3//2,totJ);continue;end
                            for J123 = 1:2:j3max
                                if !tri_check(j1,J3//2,J123//2);continue;end
                                t6j = wigner6j(Float64,j1,j2,totJ,J3p//2,J3//2,J123//2)
                                nkey = get_nkey3(J3p,J3,J123)
                                d6j_lj3[key_int][nkey] = t6j
                            end
                        end
                    end                
                end    
            end
        end
        # {t(12) t'(12) T'=1: 1/2 1/2 T123(S123)}
        d6j_dot = Dict{Int64,Dict{Int64,Float64}}() 
        Tp = 1
        for T = 1:2:3
            d6j_dot[T] = Dict{Int64,Float64}()
            tdict = d6j_dot[T]
            for t = 0:1
                for tp = 0:1
                    t6j = wigner6j(Float64,t,tp,Tp,1//2,1//2,T//2)
                    tkey = get_nkey2(t,tp)
                    tdict[tkey] = t6j
                end
            end
        end
    end
    keycg = [ zeros(Float64,3) for i=1:nthreads()]
    key6j = zeros(Float64,3);key9j = zeros(Float64,8)
    @timeit to "HOB" ninej_stuff = PreCalcHOB(paramsNN,d6j_int,to;emax_calc=e3max,mode3n=true,N3max=N3max)
    return dWS3n(dcgm0,d6j_hfspin,d6j_int,d6j_dot,d6j_lj3,d9j,d9j_int,dtri,keycg,key6j,key9j,ninej_stuff)
end 

"""
    overwrite_key6!(s12,S,s45,key)

overwrite key for 6j symbols
"""
function overwrite_key6!(s12,S,s45,key)
    key[1]=s12; key[2]=S; key[3]=s45
    return nothing
end

"""
    anti_op_isospin(params,n12,l12,s12,j12,t12,n3,l3,j3,n45,l45,s45,j45,t45,n6,l6,j6,jtot,ttot,to)
      
Function to calc. matrix element of antisymmetrizer.  
Detailed derivation can be found in e.g., Eq.(3.119) of Master Thesis by Joachim Langhammer (2010), TUDarmstadt.
"""
function anti_op_isospin(params::params3b,
                         n12,l12,s12,j12,t12,n3,l3,j3,
                         n45,l45,s45,j45,t45,n6,l6,j6,jtot,ttot,
                         to)
    dWS = params.dWS
    dict_9jHOB = dWS.ninej_stuff.dict9j_HOB
    tkeys = dWS.ninej_stuff.tkeys9j
    ex = 0.0
    if (2 * n12 + l12 + 2 * n3 + l3 != 2 * n45 + l45 + 2 * n6 + l6);return 0.0;end
    s_min = max(abs(2*s12 -1), abs(2*s45 -1))
    s_max = min(2*s12+1, 2*s45+1)
    X6 = dWS.d6j_hfspin; key6 = dWS.key6j; X9 = dWS.d9j
    for stot = s_min:2:s_max #Note: stot is already doubled
        tkey9j = tkeys[threadid()]
        lambda_min = max(abs(l12-l3),abs(l45-l6),div(abs(jtot-stot),2))
        lambda_max = min(l12+l3, l45+l6, div(jtot+stot,2))
        for lambda = lambda_min:lambda_max      
            tmp  = sqrt((2*j12+1)*(j3+1)*(2*lambda+1)*(stot+1)) 
            tmp *= sqrt((2*j45+1)*(j6+1)*(2*lambda+1)*(stot+1))
            tmp *= sqrt((2*s12+1)*(2*s45+1))     
            nkey1 = Int(get_nkey3(l12,s12,j12))
            nkey2 = Int(get_nkey2(l3,j3))
            nkey3 = Int(get_nkey3(lambda,stot,jtot))
            tmp *= X9[nkey1][nkey2][nkey3]
            if tmp == 0.0; continue;end
            nkey1 = get_nkey3(l45,s45,j45)
            nkey2 = get_nkey2(l6,j6)
            nkey3 = get_nkey3(lambda,stot,jtot)
            tmp *= X9[nkey1][nkey2][nkey3]
            tmp *= X6[get_nkey3(s12,stot,s45)]
            tmp *= HObracket(n12, l12, n3, l3, n45, l45, n6, l6, lambda, 1.0/3.0,dWS,tkey9j,dict_9jHOB,to)
            ex += tmp
        end
    end
    tmp = (-1.0)^(s12 + t12 + s45 + t45) 
    tmp *= sqrt((2*t12+1)*(2*t45+1)) 
    tmp *= X6[get_nkey3(t12,ttot,t45)]
    ex *= tmp
    anti = - 2.0 * ex / 3.0
    if (n12 == n45 && l12 == l45 && s12 == s45 && j12 == j45 && 
        t12 == t45 && n3 == n6 && l3 == l6 && j3 == j6) 
        anti = (1.0 - 2.0 * ex ) / 3.0
    end 
    return anti
end

function ret_fregvec_nonlocal(params;npow=2)
    Lambda = params.Lambda3NF
    pqs = params.meshpoints.pqs
    dim = size(pqs)[1]
    vec = zeros(Float64,dim)
    for pqidx in eachindex(params.meshpoints.wpqs)
        p,q = @view pqs[pqidx,:]
        vec[pqidx] = nonlocal_regulator_3NF(p,q,npow,Lambda)
    end
    return vec
end

"""
EGM => Eq.(19) in E.Epelbaum et al., PRC 66, 064001(2002).
Navratil => Eq.(11) in P.Navratil Few Body Syst. (2007) 41:117-140
Lambda in given in MeV
"""
function nonlocal_regulator_3NF(p,q,npow,Lambda;type="Navratil")
    if type == "EGM"
        return exp( - ((p^2 + 3/4 * q^2)/((Lambda/hc)^2))^npow)
    elseif type=="Navratil"
        return exp( - ( 0.5 * (p^2 + q^2)/((Lambda/hc)^2))^2)
    else
        @error "function nonlocal_regulator_3NF: `type`=$type is not supported!"
    end
end

"""
!PRC 66, 064001(2002)の式(A4)?
! 6 E *(4 pi)^2 δJJ' δMM' δTT' δMTMT' δl0 δλ0 δl'0 δλ0 δsj δs'j' δI1/2 δI'1/2 δtt' δss' (-1)^(t+1) w6j(1/2,1/2,t,1/2,1/2,1) 
! P.Navratil Few Body Syst. (2007) 41:117-140 式(14)に関連する表式が現れる
! <τ2・τ3> = 6 (-1)^(t+t'+T+1/2) *wigner6j(t,t',1,1/2,1/2,1/2) * wigner6j(t,t',1,1/2,1/2,T) 
"""
function IsospinDot(t_ij::Int,t_kl::Int,dT::Int,dict6j)
    fac = 6.0 * (-1)^(t_ij + t_kl + (dT+1)//2 ) 
    tkey = get_nkey2(t_ij,t_kl)
    w6j1 = dict6j[1][tkey] # T=1/2
    w6j2 = dict6j[dT][tkey] 
    r = fac * w6j1 * w6j2 * hat(t_ij) * hat(t_kl)   
    return r
end

"""
prepare R_nl(r) with 2n12+l12 + 2N + L = e12 + e3 = N3 <= N3max
Note that twice of reduced mass is used in this function
"""
function prep_Rnls(params3NF)
    N3max = params3NF.N3max
    rmass = 2*Mp*Mn/(Mp+Mn)
    b = sqrt(hc^2 /(rmass* params3NF.hw))
    ps = params3NF.meshpoints.ps
    nmesh = length(ps)
    Rnl_r = Dict{Int64,Vector{Float64}}()
    Rnl_p = Dict{Int64,Vector{Float64}}()
    for E = 0:N3max
        for n = 0:div(E,2)
            l = E - 2*n
            rnl_r = zeros(Float64,nmesh)
            rnl_p = zeros(Float64,nmesh)
            for (imesh,p) in enumerate(ps)
                rnl_p[imesh] = p * single_Rnl(p,b,n,l)
                r = p
                rnl_r[imesh] = r * single_Rnl(r,1.0/b,n,l)
            end
            tkey = get_nkey3(E,n,l)
            Rnl_p[tkey] = rnl_p
            Rnl_r[tkey] = rnl_r
        end
    end
    return Rnl_r,Rnl_p
end

