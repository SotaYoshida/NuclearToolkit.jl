"""
    make_chiEFTint(;is_show=false,itnum=1,writesnt=true,nucs=[],optimizer="",MPIcomm=false,corenuc="",ref="nucl",Operators=[])

The interface function in chiEFTint.
This generates NN-potential in momentum space and then transforms it in HO basis.
The function is exported and can be simply called make_chiEFTint() in your script. 

# Optional arguments: Note that these are mainly for too specific purposes, so you do not specify these.
- `is_show::Bool` to show `TimerOutputs`
- `itnum::Int` number of iteration for LECs calibration with HFMBPT
- `writesnt::Bool`, to write out interaction file in snt (KSHELL) format. ```julia writesnt = false``` case can be usefull when you repeat large number of calculations for different LECs.
- `nucs`
- `optimizer::String` method for LECs calibration. "MCMC","LHC","BayesOpt" are available
- `MPIcomm::Bool`, to carry out LECs sampling with HF-MBPT and affine inveriant MCMC
- `Operators::Vector{String}` specifies operators you need to use in LECs calibrations
"""
function make_chiEFTint(;is_show=false,itnum=1,writesnt=true,nucs=[],optimizer="",MPIcomm=false,corenuc="",ref="nucl",Operators=[])
    to = TimerOutput()
    io=stdout
    optHFMBPT=false
    if (optimizer!="" && nucs != []) || MPIcomm
        optHFMBPT=true; writesnt=false
    end
    if MPIcomm
        @assert optimizer == "MCMC" "when using MPI for make_chiEFTint function, optimizer should be \"MCMC\""
        @assert nucs!=[] "nucs must not be empty if you set MPIcomm=true"
        if !isdir("mpilog");run(`mkdir mpilog`);end
        MPI.Init()
        myrank = MPI.Comm_rank(MPI.COMM_WORLD)
        io = open("./mpilog/log_rank"*string(myrank)*".dat","w")
    else
        io = open("logfile.dat","w")
    end
    
    @timeit to "prep." chiEFTobj,OPTobj,d9j,HOBs = construct_chiEFTobj(optHFMBPT,itnum,optimizer,MPIcomm,io,to)
    @timeit to "NNcalc" calcualte_NNpot_in_momentumspace(chiEFTobj,to)
    @timeit to "renorm." SRG(chiEFTobj,to)
    HFdata = prepHFdata(nucs,ref,["E"],corenuc) 
    if optHFMBPT #calibrate 2n3n LECs by HFMBPT
        caliblating_2n3nLECs_byHFMBPT(itnum,optimizer,MPIcomm,chiEFTobj,OPTobj,d9j,HOBs,nucs,HFdata,to,io;Operators=Operators) 
    else # write out snt/snt.bin file
        add2n3n(chiEFTobj,to)
        @timeit to "Vtrans" dicts_tbme = TMtrans(chiEFTobj,to;writesnt=writesnt)
    end

    if io != stdout; close(io);end
    if is_show; show(to, allocations = true,compact = false);println("");end
    return true
end


"""
    ChiralEFTobject

# Fields
- `params::chiEFTparams`
- `xr_fm::Vector{Float64}` momentum mesh points in fm``^{-1}`` (0 to `pmax_fm`)
- `xr::Vector{Float64}` momentum mesh points in MeV (`xr_fm` times ``\\hbar c``)
- `wr::Vector{Float64}` weights vector for Gauss-Legendre quadrature
- `dict6j::Vector{Dict{Int64, Float64}}` dictionary of Wigner-6j symbols, `dict6j[totalJ][integer_key6j]` = 
```math
\\begin{Bmatrix} 
j_a/2&  j_b/2&   J \\\\
j_c/2&  j_d/2&   J_p
\\end{Bmatrix}
```
where `integer_key` is from the `get_nkey_from_key6j` function with ``j_a,j_b,j_c,j_d,J_p``.
- `d6j_nabla::Dict{Vector{Int64}, Float64}` dict. for Talmi-Moshinsky transformation `d6j_nabla[[ja,jb,l2,l1,0]]` = 
```math
\\begin{Bmatrix} 
j_a/2&  j_b/2&     1 \\\\
  l_2&    l_1&   1/2
\\end{Bmatrix}
```
- `d6j_int::Vector{Dict{Int64, Float64}}` dict. of Wigner-6j used for HOBs in HF-MBPT/IMSRG.
- `Rnl::Array{Float64,3}` arrays for radial functions
- `xrP_fm::Vector{Float64}` "valence" counter part of `xr_fm` (usually not used)
- `xrP::Vector{Float64}` "valence" counter part of `xr` (usually not used)
- `wrP::Vector{Float64}` "valence" counter part of `wr` (usually not used)
- `RNL::Array{Float64,3}` "valence" counter part of `Rnl` (usually not used)
- `lsjs::Vector{Vector{Vector{Int64}}}`
- `llpSJ_s::Vector{Vector{Int64}}`
- `tllsj::Vector{Int64}`
- `opfs::Vector{Vector{Float64}}`
- `ts::Vector{Float64}` mesh points for angular integrals in pion-exchange terms
- `ws::Vector{Float64}` weights vectors for angular integrals in pion-exchange terms
- `Numpn::Dict{Vector{Int64},Int64}`
- `infos::Vector{Vector{Int64}}`
- `izs_ab::Vector{Vector{Vector{Int64}}}`
- `nTBME::Int64` # of two-body matrix elements (TBMEs)
- `util_2n3n::util_2n3n`
- `LECs::LECs`
- `X9::Vector{Vector{Dict{Vector{Int64}, Float64}}}` 9j for Vtrans
- `U6::Vector{Vector{Vector{Vector{Vector{Vector{Float64}}}}}}` 6j for Vtrans
- `V12mom::Vector{Matrix{Float64}}` matrix elements of NN int in momentum space for each two-body channnel
- `V12mom_2n3n::Vector{Matrix{Float64}}` matrix elements of 2n3n int in momentum space for each two-body channnel
- `pw_channels::Vector{Vector{Int64}}` list of partial wave two-body channels like {[pnrank,L,L',S,J]}
- `dict_pwch::Vector{Dict{Vector{Int64}, Int64}}`
- `arr_pwch::Vector{Vector{Vector{Vector{Vector{Int64}}}}}`
"""
struct ChiralEFTobject
    params::chiEFTparams
    xr_fm::Vector{Float64}
    xr::Vector{Float64}
    wr::Vector{Float64}
    dict6j::Vector{Dict{Int64, Float64}}
    d6j_nabla::Dict{Vector{Int64}, Float64}  
    d6j_int::Vector{Dict{Int64, Float64}}
    Rnl::Array{Float64,3}
    xrP_fm::Vector{Float64}
    xrP::Vector{Float64}
    wrP::Vector{Float64}
    RNL::Array{Float64,3}
    lsjs::Vector{Vector{Vector{Int64}}}
    llpSJ_s::Vector{Vector{Int64}}
    tllsj::Vector{Int64}
    opfs::Vector{Vector{Float64}}
    ts::Vector{Float64}
    ws::Vector{Float64}
    Numpn::Dict{Vector{Int64},Int64}
    infos::Vector{Vector{Int64}} 
    izs_ab::Vector{Vector{Vector{Int64}}}
    nTBME::Int64
    util_2n3n::util_2n3n
    LECs::LECs
    X9::Vector{Vector{Dict{Vector{Int64}, Float64}}}
    U6::Vector{Vector{Vector{Vector{Vector{Vector{Float64}}}}}}
    V12mom::Vector{Matrix{Float64}}
    V12mom_2n3n::Vector{Matrix{Float64}}
    pw_channels::Vector{Vector{Int64}}
    dict_pwch::Vector{Dict{Vector{Int64}, Int64}} 
    arr_pwch::Vector{Vector{Vector{Vector{Vector{Int64}}}}}
end
"""

It returns 
- `chiEFTobj::ChiralEFTobject` parameters and arrays to generate NN (+2n3n) potentials. See also struct `ChiralEFoObject`.
- `OPTobj` (mutable) struct for LECs calibrations. It can be `LHSobject`/`BOobject`/`MCMCobject`/`MPIMCMCobject` struct.
- `d9j::Vector{Vector{Vector{Vector{Vector{Vector{Vector{Float64}}}}}}}` array of Wigner-9j symbols used for Pandya transformation in HF-MBPT/IMSRG calculations.
- `HOBs::Dict{Int64, Dict{Int64, Float64}}` dictionary for harmonic oscillator brackets.
"""
function construct_chiEFTobj(optHFMBPT,itnum,optimizer,MPIcomm,io,to)
    # specify chiEFT parameters
    params = init_chiEFTparams(;io=io)

    #Next, prepare momentum/integral mesh, arrays, etc.
     ## Prep WignerSymbols
    dict6j,d6j_nabla,d6j_int = PreCalc6j(params.emax,!optHFMBPT)
    ## prep. momentum mesh
    xr_fm,wr = Gauss_Legendre(0.0,params.pmax_fm,params.n_mesh); xr = xr_fm .* hc
    pw_channels,dict_pwch,arr_pwch = prepare_2b_pw_states(;io=io)    
    V12mom = [ zeros(Float64,params.n_mesh,params.n_mesh) for i in eachindex(pw_channels)]
    V12mom_2n3n = [ zeros(Float64,params.n_mesh,params.n_mesh)  for i in eachindex(pw_channels)]
    ## prep. radial functions
    rmass = Mp*Mn/(Mp+Mn)
    br = sqrt(hc^2 /(rmass* params.hw))
    Rnl = Rnl_all_ab(params,lmax,br,params.n_mesh,xr_fm)
    ## prep. for valence space oparators (usually not used)
    ntmp = ifelse(params.v_chi_order>0,params.n_mesh,3)
    xrP_fm,wrP = Gauss_Legendre(0.0,params.Pmax_fm,ntmp); xrP = xrP_fm .* hc
    RNL = Rnl_all_ab(params,lcmax,br,ntmp,xrP_fm)
    ## prep. for partial-wave decompositon
    lsjs = [[[J,J,0,J],[J,J,1,J],[J+1,J+1,1,J],[J-1,J-1,1,J],[J+1,J-1,1,J],[J-1,J+1,1,J]] for J = 0:jmax]
    llpSJ_s = [[0,0,1,1],[1,1,1,0],[1,1,0,1],[1,1,1,1],[0,0,0,0],[0,2,1,1],[1,1,1,2]]
    tllsj = zeros(Int64,5)
    opfs = [ zeros(Float64,11) for i=1:5]#T,SS,C,LS,SL terms ### why 11 or 8???
    f_ss!(opfs[2]);f_c!(opfs[3])
    ## prep. Gauss point for integrals
    ts, ws = Gauss_Legendre(-1,1,96)
    ## prep. for TBMEs        
    Numpn= Dict{Vector{Int64},Int64}()
    infos,izs_ab,nTBME = make_sp_state(params,Numpn;io=io)
    println(io,"# of channels 2bstate ",length(infos)," #TBME = $nTBME")
    ## prep. integrals for 2n3n
    util_2n3n = prep_integrals_for2n3n(params,xr,ts,ws)   
    ### specify low-energy constants (LECs)
    LECs = read_LECs(params.pottype)
    ## 9j&6j symbols for 2n (2n3n) interaction
    X9,U6 = prepareX9U6(2*params.emax)   
    chiEFTobj = ChiralEFTobject(params,xr_fm,xr,wr,dict6j,d6j_nabla,d6j_int,Rnl,
                                xrP_fm,xrP,wrP,RNL,lsjs,llpSJ_s,tllsj,opfs,ts,ws,
                                Numpn,infos,izs_ab,nTBME,util_2n3n,LECs,X9,U6,
                                V12mom,V12mom_2n3n,pw_channels,dict_pwch,arr_pwch)
      
    # make Opt stuff    
    OPTobj = prepOPT(LECs,optHFMBPT,to,io;num_cand=itnum,optimizer=optimizer,MPIcomm=MPIcomm) 
    d9j = Vector{Vector{Vector{Vector{Vector{Vector{Float64}}}}}}[ ]
    HOBs = Dict{Int64, Dict{Int64, Float64}}()
    if optHFMBPT
        d9j,HOBs = PreCalcHOB(params,d6j_int,to;io=io)
    end
    return chiEFTobj,OPTobj,d9j,HOBs
end

function calcualte_NNpot_in_momentumspace(chiEFTobj,to)
    if chiEFTobj.params.calc_NN
        OPEP(chiEFTobj,to)
        LO(chiEFTobj,to) 
        if chiEFTobj.params.chi_order >= 1
            NLO(chiEFTobj,to) 
            tpe(chiEFTobj,to) 
        end                
        if chiEFTobj.params.chi_order >= 3
            N3LO(chiEFTobj,to)
        end
        if chiEFTobj.params.chi_order >= 4
            N4LO(chiEFTobj,to)
        end
    end
    return nothing
end

function updateLECs_in_chiEFTobj!(chiEFTobj::ChiralEFTobject,targetkeys,targets)
    for (k,target) in enumerate(targetkeys)
        idx = chiEFTobj.LECs.idxs[target]
        chiEFTobj.LECs.vals[idx] = chiEFTobj.LECs.dLECs[target] = targets[k] 
    end
    return nothing
end

function add2n3n(chiEFTobj,to,it=1)
    calc_vmom_3nf(chiEFTobj,it,to)
    add_V12mom!(chiEFTobj.V12mom,chiEFTobj.V12mom_2n3n)        
    return nothing
end


function add_V12mom!(V12mom,V12mom_2n3n,a=1.0)
    for i in eachindex(V12mom)
        V12mom_2n3n[i] .+= V12mom[i] * a
    end
    return nothing
end


function genLaguerre(n,alpha,x)
    if n==0
        return 1
    elseif n==1
        return -x + (alpha+1)
    else
        s = 0.0
        for i=0:n
            bfac = gamma(n+alpha+1) / ( gamma(n-i+1) * gamma(alpha+i+1))
            s += (-1)^i * x^i / factorial(i) * bfac
        end        
        return s
    end
end

"""
    Gauss_Legendre(xmin,xmax,n;eps=3.e-16) 

Calculating mesh points `x` and weights `w` for Gauss-Legendre quadrature. This returns `x, w`.
"""
function Gauss_Legendre(xmin,xmax,n;eps=3.e-16) 
    m = div(n+1,2)
    x = zeros(Float64,n); w = zeros(Float64,n)
    xm = 0.5 * (xmin+xmax); xl = 0.5 * (xmax-xmin)
    p1 = p2 = p3 = pp = z = z1 = 0.0
    for i =1:m
        z1 = 0.0
        z = cos((pi*(i-0.25))/(n+0.5))
        hit = 0        
        while true
            hit += 1
            p1 = 1.0; p2 = 0.0
            for j=1:n
                p3=p2
                p2=p1
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j
            end
            pp = n*(z*p1-p2)/(z^2 -1.0)
            z1 = z
            z -= p1/pp
            if abs(z-z1) < eps; break;end
            if hit > 100; println("warn! in Gausse_Legendre");exit();end
        end
        x[i]=xm-xl*z
        x[n+1-i]=xm+xl*z
        w[i]=2.0*xl/( (1.0-z^2) * pp^2)
        w[n+1-i]=w[i]
    end
    return x,w
end

""" 
    QL(z,J::Int64,ts,ws)

To calculate Legendre functions of second kind by Gauss-Legendre quadrature
"""
function QL(z,J::Int64,ts,ws,QLdict)
    s = 0.0
    if J==0
        s = 0.5 * log( abs((1.0+z)/(1.0-z)) )
    elseif J==1
        s = 0.5 *z * log( abs((1.0+z)/(1.0-z)) ) -1.0
    elseif J==2
        s = 0.25 * (3.0*z^2 -1.0) * log( abs((1.0+z)/(1.0-z)) ) - 1.5*z
    elseif J==3
        s = 0.25 * (5.0*z^3 -3.0*z) * log( abs((1.0+z)/(1.0-z)) ) - 2.5*z^2 + 2.0/3.0
    else
        if J >= 0
            t = get(QLdict[J+1],z,0.0) 
            if t == 0.0
                @inbounds for (i,t) in enumerate(ts)
                    s += ws[i] * (1.0-t*t)^J / (z-t)^(J+1) 
                end 
                s *= 1.0 / 2.0^(J+1)
                QLdict[J+1][z] = s
            else
                return t
            end
        end
    end
    return s
end

"""
    make_sp_state(chiEFTobj,Numpn;io=stdout)

Defining the number of single particle states from emax and two-body channels.
"""
function make_sp_state(chiEFTobj,Numpn;io=stdout)
    emax = chiEFTobj.emax
    jab_max = 4*emax + 2
    kh = Dict{Vector{Int64},Int64}()
    kn = Dict{Vector{Int64},Int64}()
    kl = Dict{Vector{Int64},Int64}()
    kj = Dict{Vector{Int64},Int64}()
    maxsps = Int((emax+1) * (emax+2) / 2)
    println(io,"# of sp states $maxsps")    
    n = 0
    for NL=0:emax
        for L=0:NL
            if (NL-L) % 2 == 1;continue;end
            Nn= div(NL-L,2)
            for IS=-1:2:1
                jd=2*L+IS
                if jd < 0; continue;end
                n=n+1
                kh[[-1,n]]=3;  kh[[1,n]]=3
                kn[[-1,n]]=Nn; kn[[1,n]]=Nn
                kl[[-1,n]]=L;  kl[[1,n]]=L
                kj[[-1,n]]=jd; kj[[1,n]]=jd
                Numpn[[-1,Nn,L,jd]]=n
                Numpn[[ 1,Nn,L,jd]]=n
            end
        end
    end 
    infos,izs_ab,nTBME = get_twobody_channels(emax,kh,kn,kl,kj,maxsps,jab_max)
    return infos,izs_ab,nTBME
end

"""
    get_twq_2b(emax,kh,kn,kl,kj,maxsps,jab_max)

returns `infos`, `izs_ab`, `nTBME`
- `infos::Vector{Vector{Int64}}` information of two-body channeles: { [``T_z``,parity,J,dim] } 
- `izs_ab::Vector{Vecrtor{Int64}}` two-body kets: { [iza,ia,izb,ib] } where `iz*` stuff is isospin (-1 or 1) and `i*` stuff is label of sps.
For example, [-1,1,1,1] corresponds to |p0s1/2 n0s1/2>.
- `nTBME::Int` # of TBMEs
"""
function get_twobody_channels(emax,kh,kn,kl,kj,maxsps,jab_max)
    infos = [ [0,0,0,0] ];deleteat!(infos,1)
    izs_ab = [ [[0,0,0,0]] ];deleteat!(izs_ab,1)
    ichan = 0
    nTBME = 0
    for izz = -2:2:2
        for ipp = -1:2:1
            for jjx = 0:2*emax+1
                if jjx  > div(jab_max,2);continue;end
                ndim = 0
                tl = [ [0,0,0,0]];deleteat!(tl,1)
                for iza = -1:2:1
                    for izb = iza:2:1
                        if iza + izb != izz;continue;end
                        for ia = 1:maxsps
                            la = kl[[iza, ia]]; ja = kj[[iza, ia]]
                            ibmin = 1
                            if iza == izb;ibmin = ia; end                            
                            for ib = ibmin:maxsps
                                lb = kl[[izb, ib]]; jb = kj[[izb, ib]]                                
                                if (-1)^(la+lb) != ipp;continue;end
                                if tri_check(ja, jb, 2*jjx) == false; continue;end
                                if iza==izb && ia==ib && jjx % 2 == 1;continue;end
                                ndim = ndim + 1
                                push!(tl,[iza,ia,izb,ib])
                            end                            
                        end
                    end
                end
                if ndim != 0
                    ichan += 1
                    #println("ichan $ichan ndim $ndim")
                    push!(infos,[izz,ipp,jjx,ndim])
                    nTBME += Int(ndim*(ndim+1)/2)
                    push!(izs_ab,tl)
                end
            end
        end
    end
    return infos,izs_ab,nTBME
end

function def_sps_snt(emax,target_nlj)
    nljsnt = [ [0,0] ]; deleteat!(nljsnt,1)    
    tzs = [-1,1]
    dict = Dict(0=>0); delete!(dict,0)
    idx = 0; hit = 0
    for pn = 1:2
        tz = tzs[pn]
        for temax = 0:emax
            for l = temax%2:2:temax
                n = div(temax - l,2)
                jmin = 2*l-1
                if jmin < 1; jmin=1;end
                for j=jmin:2:2*l+1
                    push!(nljsnt,[n,l,j,tz])
                    idx += 1
                    for (k,tmp) in enumerate(target_nlj)
                        if tmp[1]==n && tmp[2]==l && tmp[3]==j
                            dict[idx] = k + ifelse(pn==2,length(target_nlj),0)
                        end
                    end
                end
            end
        end
    end
    return nljsnt,dict
end

function freg(p,pp,n)
    if n==0;return 1.0;end
    if n==-1 # for 3D3 pw in older EM interaction
        return 0.5 *( exp( - (p/Lambchi)^4 - (pp/Lambchi)^4 )
                      + exp( - (p/Lambchi)^6 - (pp/Lambchi)^6 ))
    end
    return exp( - (p/Lambchi)^(2*n) - (pp/Lambchi)^(2*n) )
end

"""
    read_LECs(pottype)

read LECs for a specified potential type, `em500n3lo`,`em500n3lo`,`emn500n4lo`.
"""
function read_LECs(pottype)  
    if pottype =="em500n3lo"
        return dict_em500n3lo()
    elseif pottype == "emn500n3lo"
        return dict_emn500n3lo()
    elseif pottype == "emn500n4lo"
        return dict_emn500n4lo()
    else
        @error "unknown potype $pottype"
        exit()
    end
end

"""
    calc_Vmom!(pnrank,V12mom,tdict,xr,LEC,LEC2,l,lp,S,J,pfunc,n_reg,to;is_3nf=false)    

calc. nn-potential for momentum mesh points
"""
function calc_Vmom!(chiEFTparams::chiEFTparams,pnrank,V12mom,tdict,xr,LEC,LEC2,l,lp,S,J,pfunc,n_reg,to;is_3nf=false)
    n_mesh = chiEFTparams.n_mesh
    itt = itts[pnrank]; MN = Ms[pnrank]; dwn = 1.0/MN
    V12idx = get(tdict,[itt,lp,l,S,J],-1)
    if V12idx == -1;return nothing;end
    V = V12mom[V12idx]
    @inbounds for i= 1:n_mesh
        x = xr[i]; ex = sqrt(1.0+(dwn.*x)^2)
        @inbounds for j = 1:n_mesh
            y = xr[j]; ey = sqrt(1.0+(dwn.*y)^2)
            ree = 1.0/sqrt(ex*ey)
            if is_3nf;ree=1.0;end
            fac = pfunc(x,y,LEC,LEC2) * freg(x,y,n_reg) *ree
            V[i,j] += fac
        end
    end
    return nothing
end

"""
    Rnl_all_ab(lmax,br,n_mesh,xr_fm)

Returns array for radial functions (prop to generalized Laguerre polynomials) HO w.f. in momentum space.
Rnlk(l,n,k)=sqrt(br) * R(n,L,Z) *Z with Z=br*k (k=momentum in fm^-1)
"""
function Rnl_all_ab(chiEFTobj,lmax_in,br,n_mesh,xr_fm)
    Nnmax = chiEFTobj.Nnmax
    Rnl = zeros(Float64,Nnmax+1,lmax_in+1,n_mesh)
    for l=0:lmax_in
        for kidx=1:n_mesh
            pb = br * xr_fm[kidx]
            pb2 = pb^2
            fexp = exp(-0.5*pb2)
            fpow = pb^(l+1)
            for n=0:Nnmax
                fac = sqrt(2.0*br*factorial(n) / gamma(n+l+3/2)) * fexp * fpow
                fac *= genLaguerre(n,l+1//2,pb2)
                Rnl[n+1,l+1,kidx]= fac
            end
        end
    end
    return Rnl
end

"""
    prepare_2b_pw_states(;io=stdout)

preparing two-body channels in terms of <Lp,S,J| |L,S,J>.
For example, [pnrank,L,Lp,S,J] = [0, 0, 2, 1, 1] corresponds to proton-neutron 3S1-3D1 channel.
"""
function prepare_2b_pw_states(;io=stdout)
    #iz,lz1,lz2,isz,jz
    pw_channels = [ [0,0,0,0,0] ];deleteat!(pw_channels,1)
    dict_pwch = [Dict{Vector{Int64},Int64}() for pnrank=1:3] 
    arr_pwch = [[[[ zeros(Int64,j+iss-abs(j-iss)+1) for ll1=abs(j-iss):j+iss ] for j=0:jmax ] for iss=0:1 ] for pnrank=1:3]    
    num=0
    ## pp iz = -2
    itt = 1
    for iss=0:1
        for j=0:jmax
            for ll1=abs(j-iss):j+iss
                for ll2=abs(j-iss):j+iss
                    if (-1)^ll1 != (-1)^ll2;continue;end
                    if itt != Int( (1+(-1)^(ll1+iss))/2);continue;end
                    if (-1)^(ll1+iss+itt) != -1;continue;end
                    num=num+1
                    push!(pw_channels,[-2,ll1,ll2,iss,j])
                    dict_pwch[1][[-2,ll1,ll2,iss,j]] = num
                    arr_pwch[1][iss+1][j+1][ll1-abs(j-iss)+1][ll2-abs(j-iss)+1] = num
                end
            end
        end
    end
    ##  nn iz = 2
    pnrank=3
    for iss=0:1
        for j=0:jmax
            for ll1=abs(j-iss):j+iss
                for ll2=abs(j-iss):j+iss
                    if (-1)^ll1 != (-1)^ll2;continue;end
                    if itt != Int( (1+(-1)^(ll1+iss))/2);continue;end
                    if (-1)^(ll1+iss+itt) != -1;continue;end
                    num=num+1
                    push!(pw_channels,[2,ll1,ll2,iss,j])
                    dict_pwch[pnrank][[2,ll1,ll2,iss,j]] = num
                    arr_pwch[pnrank][iss+1][j+1][ll1-abs(j-iss)+1][ll2-abs(j-iss)+1] = num
                end
            end
        end
    end
    ## pn iz = 0
    pnrank=2
    for iss=0:1
        for ktt=1:2
            if (iss==0 && ktt==1) || (iss==1 && ktt==2);itt=1;end
            if (iss==0 && ktt==2) || (iss==1 && ktt==1);itt=0;end
            for j=0:jmax
                for ll1=abs(j-iss):j+iss
                    for ll2=abs(j-iss):j+iss
                        if (-1)^ll1 != (-1)^ll2;continue;end
                        if itt != Int((1+(-1)^(ll1+iss))/2) ;continue;end
                        if (-1)^(ll1+iss+itt) != -1;continue;end
                        num=num+1
                        push!(pw_channels,[0,ll1,ll2,iss,j])
                        dict_pwch[pnrank][[0,ll1,ll2,iss,j]] = num
                        arr_pwch[pnrank][iss+1][j+1][ll1-abs(j-iss)+1][ll2-abs(j-iss)+1] = num
                    end
                end
            end
        end
    end
    println(io,"# of two-body states $num")
    return pw_channels,dict_pwch,arr_pwch
end

function calc_coulomb(chiEFTobj,xs,ws,Vcoulomb,pw_channels,num2bch,Rnl;meshp=100)
    hw = chiEFTobj.hw
    Nnmax = chiEFTobj.Nnmax
    rmass = 0.5 * Mp
    brange = sqrt(hc2 / (rmass*hw))
    rmax = 20.0
    r,w = Gauss_Legendre(0.0,rmax,meshp)
    ra1 = zeros(Float64,meshp);rb1 = zeros(Float64,meshp)
    ra2 = zeros(Float64,meshp);rb2 = zeros(Float64,meshp)
    rnl1 = zeros(Float64,meshp)
    rnl2 = zeros(Float64,meshp)
    vcl = zeros(Float64,meshp)
    fmcoul(vcl,r)
    memo1=[0]; memo2=[0]
    for num=1:num2bch
        vcoul = Vcoulomb[num]
        iz12,l1,l2,isz,jj = pw_channels[num]
        if iz12 !=-2; continue;end
        if l1 !=l2; continue;end
        if jj > jmax;continue;end
        nmax = Nnmax
        memo1[1]=0              
        for n1=0:nmax
            howf(num,brange,n1,l1,r,ra1,rb1,rnl1,meshp,memo1) # => rnl1
            memo2[1]=0
            for n2=0:n1
                howf(num,brange,n2,l2,r,ra2,rb2,rnl2,meshp,memo2) # => rnl2
                scoul = 0.0
                for i=1:meshp
                    scoul += rnl1[i]*rnl2[i]*w[i]*vcl[i]
                end
                vcoul[n1+1,n2+1]= vcoul[n2+1,n1+1]=scoul
            end
        end
    end
    return nothing
end
function fmcoul(vcl,r)
    for (i,x) in enumerate(r)
        vcl[i] = hc *fsalpha /x
    end
    return nothing
end

function howf(num,brange,n,l,r,ra,rb,ret_rnl,meshp,memo)
    #BR,N,L,R,RNL,NMESH
    kmax = l +1
    gam = sqrt(pi)
    for k=1:kmax
        gam *= (0.5+k-1)
    end
    if memo[1] < 1
        for i = 1:meshp
            z = r[i]/brange
            zz = z*z
            zl = 1.5 + l
            ra[i] = sqrt(2.0/(gam*brange)) * z^(l+1) *exp(-0.5*zz)
            rb[i] = ra[i] * (zl-zz)/sqrt(zl)
        end
        memo[1] = 1
    end
    if n < 1; ret_rnl .= ra; return nothing;end
    if n ==1; ret_rnl .= rb; return nothing;end

    if n == memo[1]-1
        ret_rnl .= ra
        return nothing
    end
    if n == memo[1]
        ret_rnl .= rb
        return nothing
    end
    if n > memo[1]
        memo[1] += 1
        a = memo[1]*2+l -0.5
        b = sqrt(memo[1] *(memo[1]+l+0.5))
        c = sqrt( (memo[1]-1)*(memo[1]+l-0.5))
        for i =1:meshp
            zz = (r[i]/brange)^2
            ret_rnl[i] = ((a-zz)*rb[i]-c*ra[i])/b
            ra[i] = rb[i]
            rb[i] = ret_rnl[i]
        end
        while memo[1] < n
            memo[1] += 1
            a = memo[1]*2+l -0.5
            b = sqrt(memo[1] *(memo[1]+0.5))
            c = sqrt( (memo[1]-1)*(memo[1]+l-0.5))
            for i =1:meshp
                zz = (r[i]/brange)^2
                ret_rnl[i] = ((a-zz)*rb[i]-c*ra[i])/b
                ra[i] = rb[i]
                rb[i] = ret_rnl[i]
            end
        end
        return nothing
    end
    if n < memo[1] -1 
        memo[1] = 1
    end
    return nothing
end


function Vrel(chiEFTobj,to) 
    V12mom = chiEFTobj.V12mom; pw_channels = chiEFTobj.pw_channels; Rnl = chiEFTobj.Rnl
    xr_fm = chiEFTobj.xr_fm; wr = chiEFTobj.wr; n_mesh = chiEFTobj.params.n_mesh
    Nnmax= chiEFTobj.params.Nnmax; num2bch = length(pw_channels)
    V12ab = [zeros(Float64,Nnmax+1,Nnmax+1) for i=1:num2bch]
    Vcoulomb = [zeros(Float64,Nnmax+1,Nnmax+1) for i=1:num2bch]
    if chiEFTobj.params.coulomb
        calc_coulomb(chiEFTobj.params,xr_fm,wr,Vcoulomb,pw_channels,num2bch,Rnl)
    end
    x = zeros(Float64,Nnmax+1,n_mesh)
    for num = 1:num2bch
        Vtmp = V12mom[num]
        Vab = V12ab[num]
        vcoul = Vcoulomb[num]
        iz,l1,l2,isz,jz = pw_channels[num]
        @inbounds for n1 = 0:Nnmax
            tR = @views Rnl[n1+1,l1+1,:]
            tx = @views x[n1+1,:]
            @inbounds for k2=1:n_mesh
                sum=0.0
                @inbounds for k1=1:n_mesh
                    vk1k2=Vtmp[k1,k2]
                    if vk1k2 == 0.0; continue;end
                    pkx1=tR[k1]
                    sum += pkx1*vk1k2*xr_fm[k1]*wr[k1]
                end
                tx[k2]=sum*wr[k2]*xr_fm[k2]
            end
        end
        @inbounds for n1=0:Nnmax
            @inbounds for n2=0:Nnmax
                tR = @views Rnl[n2+1,l2+1,:]
                vsum=0.0
                @inbounds for k2=1:n_mesh
                    pkx2=tR[k2]
                    vsum += x[n1+1,k2]*pkx2
                end
                phase=(-1.0)^(n1+n2) 
                t_vcoul = vcoul[n1+1,n2+1]
                Vab[n1+1,n2+1]= phase * vsum + t_vcoul
            end
        end
    end
    return V12ab
end


"""
    hw_formula(A,fnum)

empirical formula for harmonis oscillator parameter hw by mass number A
fnum=2: for sd-shell, Ref. J. Blomqvist and A. Molinari, Nucl. Phys. A106, 545 (1968).
"""
function hw_formula(A,fnum)
    hw = 0.0
    if fnum == 2        
        hw = 45.0 * (A^(-1.0/3.0)) -25.0 * (A^(-2.0/3.0))
    else
        hw = 41.0 * (A^(-1.0/3.0))
    end
    return hw
end
