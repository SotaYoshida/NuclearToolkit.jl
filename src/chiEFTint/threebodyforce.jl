struct dWS3N
    dtri::Dict{Int64,Float64}
    dcgm0::Dict{Int64,Float64}
    d6j_int::Dict{UInt64,Float64}
    d6j_lj::Dict{UInt64,Float64}
    d9j_int::Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}
    d9j_lsj::Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}
    dictHOB::Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}
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
    e2max::Int64
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
    dWS::dWS3N    
    main_basis::String
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

function dev_param3NF()
    TFrenorm = false
    main_basis="Jacobi"   #K.Hebeler    
    main_basis="JacobiHO" #P.Navratil  
    regulator = "nonlocal"
    regulator = "local"
    Np = Nq = 10 #25
    Nangle = 10 #25
    pmax3 = 10.0
    rmax3 = 10.0
    e3max = 4
    N3max = 10 # to be 40?
    j3max = 9 # said to be sufficient, but can be e3max*2 + 3

    # for check
    N3max = 5
    j3max = 9
    e3max = 3

    L3max = min(N3max,div(j3max-1,2))
    J12max = 2*j3max
    chiorder3b = 2 # NNLO
    Lambda3NF = 400
    LambdaChi = 700.0
    return Np,Nq,Nangle,e3max,N3max,L3max,j3max,J12max,hw,chiorder3b,Lambda3NF,LambdaChi,TFrenorm,pmax3,rmax3,main_basis,regulator
end

function get_ket_JacobiHO(n,l,s,j,t,N,L,J,dJ3,dT3)
    e1 = 2*n+l
    e2 = 2*N+L
    Eket = e1+e2
    return ket_JacobiHO(n,l,s,j,t,N,L,J,dJ3,dT3,e1,e2,Eket)
end

function test3NF(;param_str="dev",target_LECs=["c1_NNLO","c3_NNLO","c4_NNLO","cD","cE"],fn_params="optional_parameters.jl",
                is_show=false)
    to = TimerOutput()
    io =  select_io(false,"",[])
    paramsNN = init_chiEFTparams(;io=io,fn_params=fn_params)
    LECs = read_LECs(paramsNN.pottype)
    params3N = get_params3N(param_str,paramsNN,LECs,to)
    if params3N.main_basis == "Jacobi" # K.Hebeler type                
        Calc_3NF_in_Jacobi_coordinate(params3N,LECs,to)
    elseif params3N.main_basis == "JacobiHO" # P.Navratil type
        Calc_3NF_in_JacobiHO_coorrdinate(params3N,LECs,to)
    else
        @error "main_basis=$main_basis is not supported"
    end
    if is_show 
        show(to, allocations = true,compact = false);println("")
    end
end

function make_3b_mesh_mats(pmax3,Np,Nq,Nangle)
    ps,wps = Gauss_Legendre(0.0,pmax3,Np)
    qs,wqs = Gauss_Legendre(0.0,pmax3,Nq)
    pqs,wpqs,dim = make_pq(ps,wps,qs,wqs)
    us,wus = Gauss_Legendre(-1.0,1.0,Nangle)
    return mesh_3NF(ps,wps,qs,wqs,pqs,wpqs,us,wus)
end

function get_params3N(param_str,paramsNN,LECs,to;io=stdout)    
    e2max = paramsNN.emax
    Np,Nq,Nangle,e3max,N3max,L3max,j3max,J12max,hw,chiorder3b,Lambda3NF,LambdaChi,TFrenorm,pmax3,rmax3,main_basis,regulator = dev_param3NF()
    if param_str != "dev"
        @error "param_str = $paramstr is not supported now"
    end

    dWS = prep_dWS3N(N3max,J12max,j3max,to)
    
    meshpoints = make_3b_mesh_mats(pmax3,Np,Nq,Nangle)
    return params3b(e2max,e3max,N3max,L3max,j3max,J12max,hw,chiorder3b,Lambda3NF,LambdaChi,TFrenorm,pmax3,rmax3,meshpoints,LECs,dWS,main_basis,regulator)
end

function overwrite3NFparams(fn,e3max_in,N3max_in,L3max_in,j3max_in,J12max_in,hw_in,chiorder3b_in,Lambda3NF_in,LambdaChi_in,TFrenorm_in,pmax3_in,rmax3_in)
    e3max_r = e3max_in
    if @isdefined(e3max); e3max_r = e3max; end
    N3max_r = N3max_in
    if @isdefined(e3max); _r = e3max; end

    L3max_r = L3max_in
    j3max_r = j3max_in
    J12max_r = J12max_in
    hw_r = hw_in
    chiorder3b_r = chiorder3b_in
    Lambda3NF_r = Lambda3NF_in
    LambdaChi_r = LambdaChi_in
    TFrenorm_r = TFrenorm_in
    pmax3_r = pmax3_in
    rmax3_r = rmax3_in

    return e3max_r,N3max_r,L3max_r,j3max_r,J12max_r,hw_r,chiorder3b_r,Lambda3NF_r,LambdaChi_r,TFrenorm_r,pmax3_r,rmax3_r
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
    if (2 * n12 + l12 + 2 * n3 + l3 != 2 * n45 + l45 + 2 * n6 + l6);return 0.0;end

    ex = 0.0
    d6j = params.dWS.d6j_lj
    d9j = params.dWS.d9j_lsj
    dictHOB = params.dWS.dictHOB

    s_min = max(abs(2*s12 -1), abs(2*s45 -1))
    s_max = min(2*s12+1, 2*s45+1)
    hit = 0
    for stot = s_min:2:s_max #Note: stot is already doubled
        lambda_min = max(abs(l12-l3),abs(l45-l6),div(abs(jtot-stot),2))
        lambda_max = min(l12+l3, l45+l6, div(jtot+stot,2))
        for lambda = lambda_min:lambda_max
            tmp  = sqrt((2*j12+1)*(j3+1)*(2*lambda+1)*(stot+1)) 
            tmp *= sqrt((2*j45+1)*(j6+1)*(2*lambda+1)*(stot+1))
            tmp *= sqrt((2*s12+1)*(2*s45+1))     
            tmp *= call_d9j_lsj(l12*2,s12*2,j12*2,l3*2,1,j3,2*lambda,stot,jtot,d9j)     
            if tmp == 0.0; continue;end
            tmp *= call_d9j_lsj(l45*2,s45*2,j45*2,l6*2,1,j6,2*lambda,stot,jtot,d9j)        
            t6j = call_d6j(1,1,s12*2,1,stot,s45*2,d6j)
            tmp *= t6j
            tmp *= get_dictHOB(n12, l12, n3, l3, n45, l45, n6, l6, lambda, dictHOB)
            ex += tmp
            hit += 1
        end
    end
    tmp = (-1.0)^(s12 + t12 + s45 + t45) 
    tmp *= sqrt((2*t12+1)*(2*t45+1)) 
    tmp *= call_d6j(1,1,t12*2,1,ttot,t45*2,d6j)
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
<τ2・τ3> = 6 (-1)^(t+t'+T+1/2) *wigner6j(t,t',1,1/2,1/2,1/2) * wigner6j(t,t',1,1/2,1/2,T) 
"""
function IsospinDot(t_ij::Int,t_kl::Int,dT::Int,d6j)
    fac = 6.0 * (-1)^(t_ij + t_kl + div(dT+1,2) ) 
    if t_ij == t_kl == 0; return 0.0;end
    w6j1 = call_d6j(t_ij*2,t_kl*2,2,1,1,1,d6j)
    w6j2 = call_d6j(t_ij*2,t_kl*2,2,1,1,dT,d6j)
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

