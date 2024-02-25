"""
    make_chiEFTint(;is_show=false,itnum=1,writesnt=true,nucs=[],optimizer="",MPIcomm=false,corenuc="",ref="nucl",Operators=[],fn_params="optional_parameters.jl",write_vmom=false,do_svd=false)

The interface function in chiEFTint.
This generates NN-potential in momentum space and then transforms it in HO basis to give inputs for many-body calculations.
The function is exported and can be simply called make_chiEFTint() in your script. 

# Optional arguments: Note that these are mainly for too specific purposes, so you do not specify these.
- `is_show::Bool` to show `TimerOutputs`
- `itnum::Int` number of iteration for LECs calibration with HFMBPT
- `writesnt::Bool`, to write out interaction file in snt (KSHELL) format. ```julia writesnt = false``` case can be usefull when you iteratively calculate with different LECs.
- `nucs` target nuclei used for LECs calibration with HFMBPT
- `optimizer::String` method for LECs calibration. "MCMC","LHS","BayesOpt" are available
- `MPIcomm::Bool`, to carry out LECs sampling with HF-MBPT and affine inveriant MCMC
- `Operators::Vector{String}` specifies operators you need to use in LECs calibrations
- `fn_params::String` path to file specifying the optional parameters
- `write_vmom::Bool` to write out in vmom partial wave channels
"""
function make_chiEFTint(;is_show=false,itnum=1,writesnt=true,nucs=[],optimizer="",MPIcomm=false,corenuc="",ref="nucl",Operators=[],fn_params="optional_parameters.jl",write_vmom=false,do_svd=false,do2n3ncalib=false)
    to = TimerOutput()    
    if (optimizer!="" && nucs != []) || MPIcomm; do2n3ncalib=true; writesnt=false; end
    io = select_io(MPIcomm,optimizer,nucs)
    @timeit to "prep." chiEFTobj,OPTobj,dWS = construct_chiEFTobj(do2n3ncalib,itnum,optimizer,MPIcomm,io,to;fn_params)
    @timeit to "NNcalc" calcualte_NNpot_in_momentumspace(chiEFTobj,to)
    BE_d_bare = Calc_Deuteron(chiEFTobj,to;io=io)
    @timeit to "renorm." SRG(chiEFTobj,to)
    BE_d_srg = Calc_Deuteron(chiEFTobj,to;io=io)
    if chiEFTobj.params.srg
        println("E(2H): bare = ",@sprintf("%9.6f", BE_d_bare),
                " srg = ", @sprintf("%9.6f", BE_d_srg), " Diff.", @sprintf("%8.3e", BE_d_bare - BE_d_srg))
    else
        println("E(2H): bare = ",@sprintf("%9.6f", BE_d_bare))        
    end
    HFdata = prepHFdata(nucs,ref,["E"],corenuc)

    if do_svd; target_LSJ = [[0,0,0,0],[0,2,1,1],[1,1,1,0],[2,2,0,2]]; svd_vmom(chiEFTobj,target_LSJ); end

    if write_vmom
        target_LSJ = [[0,0,1,1],[1,1,1,0],[1,1,0,1],[1,1,1,1],[0,0,0,0],[0,2,1,1],[3,3,1,3]]
        write_onshell_vmom(chiEFTobj,2,target_LSJ;label="pn"); write_onshell_vmom(chiEFTobj,3,target_LSJ;label="nn")
        #momplot(chiEFTobj,2,target_LSJ; fnlabel=ifelse(chiEFTobj.params.srg,"srg","bare"))
        #momplot(chiEFTobj,3,target_LSJ; fnlabel=ifelse(chiEFTobj.params.srg,"srg","bare"))
    end

    if do2n3ncalib #calibrate 2n3n LECs by HFMBPT
        caliblating_2n3nLECs_byHFMBPT(itnum,optimizer,MPIcomm,chiEFTobj,OPTobj,dWS,nucs,HFdata,to,io;Operators=Operators)        
    else # write out snt/snt.bin file
        calc_vmom_3nf(chiEFTobj,1,to)
        if chiEFTobj.params.calc_EperA; calc_nuclearmatter_in_momspace(chiEFTobj,to,io);end
        @timeit to "Vtrans" dicts_tbme = TMtrans(chiEFTobj,dWS,to;writesnt=writesnt)
    end
    if io != stdout; close(io);end
    show_TimerOutput_results(to;tf=is_show)
    
    return true
end

"""
    construct_chiEFTobj(do2n3ncalib,itnum,optimizer,MPIcomm,io,to;fn_params="optional_parameters.jl")

It returns 
- `chiEFTobj::ChiralEFTobject` parameters and arrays to generate NN (+2n3n) potentials. See also struct `ChiralEFoObject`.
- `OPTobj` (mutable) struct for LECs calibrations. It can be `LHSobject`/`BOobject`/`MCMCobject`/`MPIMCMCobject` struct.
"""
function construct_chiEFTobj(do2n3ncalib,itnum,optimizer,MPIcomm,io,to;fn_params="optional_parameters.jl")
    # specify chiEFT parameters
    params = init_chiEFTparams(;io=io,fn_params=fn_params)
    @timeit to "prep dWS" dWS = prep_dWS2n(params,to)

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
    tllsj = zeros(Int64,5)
    opfs = [ zeros(Float64,11) for i=1:5]#T,SS,C,LS,SL terms 
    f_ss!(opfs[2]);f_c!(opfs[3])
    ## prep. Gauss point for integrals
    ts, ws = Gauss_Legendre(-1,1,96)
    ## prep. for TBMEs        
    infos,izs_ab,nTBME = make_sp_state(params;io=io)
    println(io,"# of channels 2bstate ",length(infos)," #TBME = $nTBME")
    ## prep. integrals for 2n3n
    @timeit to "util2n3n" util_2n3n = prep_integrals_for2n3n(params,xr,ts,ws,to)
    ### specify low-energy constants (LECs)
    LECs = read_LECs(params.pottype)
    chiEFTobj = ChiralEFTobject(params,xr_fm,xr,wr,Rnl,
                            xrP_fm,xrP,wrP,RNL,lsjs,tllsj,opfs,ts,ws,
                            infos,izs_ab,nTBME,util_2n3n,LECs,
                            V12mom,V12mom_2n3n,pw_channels,dict_pwch,arr_pwch)
    # make Opt stuff    
    @timeit to "OPTobj" OPTobj = prepOPT(LECs,do2n3ncalib,to,io;num_cand=itnum,optimizer=optimizer,MPIcomm=MPIcomm) 
    return chiEFTobj,OPTobj,dWS
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

function updateLECs_in_chiEFTobj!(chiEFTobj::ChiralEFTobject,targetkeys,targets::Vector{Float64})
    for (k,target) in enumerate(targetkeys)
        idx = chiEFTobj.LECs.idxs[target]
        chiEFTobj.LECs.vals[idx] = chiEFTobj.LECs.dLECs[target] = targets[k] 
    end
    return nothing
end

function add2n3n(chiEFTobj::ChiralEFTobject,to,it=1)
    calc_vmom_3nf(chiEFTobj,1,to)
    add_V12mom!(chiEFTobj.V12mom,chiEFTobj.V12mom_2n3n)        
    return nothing
end

function add_V12mom!(V12mom::VM,V12mom_2n3n::VM,a=1.0) where VM<:Vector{Matrix{Float64}}
    for i in eachindex(V12mom)
        V12mom_2n3n[i] .+= V12mom[i] * a
    end
    return nothing
end
"""
    genLaguerre(n::Int,alpha,x)
returns generalized Laguaerre polynomials, ``L^{\\alpha}_n(x)``
"""
function genLaguerre(n::Int,alpha,x)
    if n==0
        return 1.0
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
    Legendre(n,x)
function to calculate Legendre polynomial ``P_n(x)``
"""
function Legendre(n,x)
    nmax = Int(ifelse(n%2==0,n//2,(n-1)//2))
    tsum = 0.0
    for k = 0:nmax
        tsum += (-1)^k * factorial(big(2*n-2*k)) / factorial(big(k)) / factorial(big(n-k)) / factorial(big(n-2*k)) *x^(n-2*k)
    end
    ret = tsum / 2^n
    return Float64(ret)
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
            if hit > 100; println("warn! in Gauss_Legendre");exit();end
        end
        x[i]=xm-xl*z
        x[n+1-i]=xm+xl*z
        w[i]=2.0*xl/( (1.0-z^2) * pp^2)
        w[n+1-i]=w[i]
    end
    @assert abs(sum(w) - (xmax-xmin)) < 1.e-8 "GL integral $(sum(w)) not match xmax-xmin $(xmax-xmin)"
    return x,w
end

""" 
    QL(z,J::Int64,ts,ws)

To calculate Legendre functions of second kind, which are needed for pion-exchange contributions, by Gauss-Legendre quadrature.
"""
function QL(z,J::Int64,ts,ws,QLdict;zthreshold=1.0)
    if J < 0; return 0.0;end
    val = 0.0
    if z > zthreshold 
        val =  QL_numeric_dict(z,J,ts,ws,QLdict)
    else
        println("z $z J $J recursive read")
        val = QL_recursive(z,J)
    end
    return val
end
function QL_numeric_dict(z,J,ts,ws,QLdict)
    ret = get(QLdict[J+1],z,0.0)
    if ret == 0.0
        ret = QL_numeric(z,J,ts,ws)
        QLdict[J+1][z] = ret
    end
    return ret
end
function QL_numeric(z,J,ts,ws)
    ret = 0.0
    @inbounds for (i,t) in enumerate(ts)
        ret += ws[i] * (1.0-t*t)^J / (2*(z-t))^(J+1)
    end 
    return ret
end
function QL_numeric_fac(z,J,ts,ws,factor_vec::Vector{Float64})
    ret = 0.0
    @inbounds for (i,t) in enumerate(ts)
        ret += factor_vec[i] * ws[i] * (1.0-t*t)^J / (2*(z-t))^(J+1) 
    end 
    return ret
end
function QL_recursive(z,J)
    s = 0.0
    logpart = log( abs((1.0+z)/(1.0-z)) )
    Q0 = 0.5 * logpart
    Q1 = 0.5 * z * logpart -1.0
    if J==0
        s = Q0
    elseif J==1
        s = Q1 
    else
        QJ = Q1; QJm = Q0; tJ = 1
        while tJ < J
            tQ = ((2*tJ+1) * z * QJ - tJ * QJm ) /(tJ+1)
            QJm = QJ; QJ = tQ
            tJ += 1
        end
        s = QJ
    end
    return s
end

"""
    make_sp_state(chiEFTparams;io=stdout)

Defining the two-body channels or kets.
"""
function make_sp_state(chiEFTobj::chiEFTparams;io=stdout)
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
                n += 1
                kh[[-1,n]]=3;  kh[[1,n]]=3
                kn[[-1,n]]=Nn; kn[[1,n]]=Nn
                kl[[-1,n]]=L;  kl[[1,n]]=L
                kj[[-1,n]]=jd; kj[[1,n]]=jd
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
function get_twobody_channels(emax::Int,kh,kn,kl,kj,maxsps::Int,jab_max::Int)
    infos = Vector{Int64}[ ]  
    izs_ab = Vector{Vector{Int64}}[ ] 
    ichan = nTBME = 0
    num2b = 0
    for izz = -2:2:2
        for ipp = -1:2:1
            for jjx = 0:2*emax+1
                if jjx  > div(jab_max,2);continue;end
                ndim = 0
                tl = Vector{Int64}[ ]
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
                    num2b +=div(ndim*(ndim+1),2)
                    push!(infos,[izz,ipp,jjx,ndim])
                    nTBME += Int(ndim*(ndim+1)/2)
                    push!(izs_ab,tl)
                end
            end
        end
    end
    # num1b = 0
    # nljs = Vector{Int64}[ ]
    # for e = 0:emax
    #     for n = 0:div(e,2)
    #         l = e - 2*n
    #         for j = max(1,2*l-1):2:2*l+1
    #             push!(nljs,[n,l,j])
    #         end
    #     end
    # end
    # for ia = 1:length(nljs)
    #     oa = nljs[ia]
    #     for ib = ia+1:length(nljs)
    #         ob = nljs[ib]
    #         if oa[2]==ob[2] && oa[3]==ob[3]
    #             num1b += 1
    #         end
    #     end
    # end
    # numtot = 2* num1b + num2b
    # println("num tot $numtot num1b $num1b num2b $num2b")
    return infos,izs_ab,nTBME
end

"""
    def_sps_snt(emax,target_nlj)
Defining dicts for single particle states.
One can truncate sps by specifying `target_nlj`, but it is usually empty so this function is usually not used.
"""
function def_sps_snt(emax,target_nlj)
    nljsnt = Vector{Int64}[ ] 
    tzs = [-1,1]
    dict = Dict{Int64,Int64}()
    idx = 0
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

"""
    freg(p,pp,n;Lambchi=500.0)
the regulator function used for NN or 2n3n contribution, Eq.(4.63) in EM's review:
``f(p,p') = \\exp [ -(p/\\Lambda)^{2n} -(p'/\\Lambda)^{2n}  ]``
"""
function freg(p,pp,n;Lambchi=500.0)
    if n==0;return 1.0;end
    if n==-1 # for 3D3 pw in older EM interaction
        return 0.5 *( exp( - (p/Lambchi)^4 - (pp/Lambchi)^4 )
                      + exp( - (p/Lambchi)^6 - (pp/Lambchi)^6 ))
    end
    return exp( - (p/Lambchi)^(2*n) - (pp/Lambchi)^(2*n) )
end

"""
    read_LECs(pottype)

read LECs for a specified potential type, `em500n3lo`,`em500n3lo`,`emn500n4lo`, `nnlosat`.
"""
function read_LECs(pottype)  
    if pottype =="em500n3lo"
        return dict_em500n3lo()
    elseif pottype == "emn500n3lo"
        return dict_emn500n3lo()
    elseif pottype == "emn500n4lo"
        return dict_emn500n4lo()
    elseif pottype == "nnlosat" 
        return dict_nnlosat()
    else
        @error "unknown potype $pottype"       
    end
end

"""
    calc_Vmom!(pnrank,V12mom,tdict,xr,LEC,LEC2,l,lp,S,J,pfunc,to;is_3nf=false)    

calc. NN-potential for momentum mesh points
"""
function calc_Vmom!(params::chiEFTparams,pnrank,V12mom,tdict,xr,LEC,LEC2,l,lp,S,J,pfunc,n_reg,to;is_3nf=false)
    n_mesh = params.n_mesh
    itt = 2 *(pnrank -2); MN = Ms[pnrank]; dwn = 1.0/MN
    V12idx = get(tdict,[itt,lp,l,S,J],-1)
    if V12idx == -1;return nothing;end
    V = V12mom[V12idx]
    @inbounds for i= 1:n_mesh
        x = xr[i]; ex = sqrt(1.0+(dwn.*x)^2)
        @inbounds for j = 1:n_mesh
            y = xr[j]; ey = sqrt(1.0+(dwn.*y)^2)
            ree = 1.0/sqrt(ex*ey)
            if is_3nf;ree=1.0;end
            fac = pfunc(x,y,LEC,LEC2) * freg(x,y,n_reg;Lambchi=params.Lambda_cutoff) *ree
            V[i,j] += fac
        end
    end
    return nothing
end

"""
    Rnl_all_ab(lmax,br,n_mesh,xr_fm)

Returns array for radial functions (prop to generalized Laguerre polynomials) HO w.f. in momentum space.
`Rnlk(l,n,k)=sqrt(br) * R(n,L,Z) *Z` with `Z=br*k` (k is momentum in fm``{}^{-1}``)
"""
function Rnl_all_ab(chiEFTobj,lmax_in,br,n_mesh,xr_fm)
    Nnmax = chiEFTobj.Nnmax
    Rnl = zeros(Float64,Nnmax+1,lmax_in+1,n_mesh)
    for l=0:lmax_in
        for kidx=1:n_mesh
            pb = br * xr_fm[kidx]; pb2 = pb^2
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
simply calculate ``R_{nl}(p,b) = \\sqrt{ \\frac{2 n! b^3}{\\Gamma(n+l+3/2)} (pb)^l e^{-p^2b^2/2} L^{l+1/2}_{n}(p^2b^2) }``
"""
function single_Rnl(p,b,n,l)
    pb  = p*b
    return sqrt(2.0*factorial(n) * b^3 / gamma(n+l+3/2)) * (pb^l) *exp(-0.5*pb^2) * genLaguerre(n,l+1//2,pb^2)
end

"""
    prepare_2b_pw_states(;io=stdout)

preparing two-body partical-wave channels, <Lp,S,J| |L,S,J>.
For example, [pnrank,L,Lp,S,J] = [0, 0, 2, 1, 1] corresponds to proton-neutron 3S1-3D1 channel.
"""
function prepare_2b_pw_states(;io=stdout)
    pw_channels = Vector{Int64}[ ]
    dict_pwch = [Dict{Vector{Int64},Int64}() for pnrank=1:3] 
    arr_pwch = [[[[ zeros(Int64,j+iss-abs(j-iss)+1) for ll1=abs(j-iss):j+iss ] for j=0:jmax ] for iss=0:1 ] for pnrank=1:3]    
    num=0
    for pnrank = [1,3,2]
        itt = 1
        Tz = 2*pnrank -4 
        for S=0:1
            for j=0:jmax
                for ll1=abs(j-S):j+S
                    for ll2=abs(j-S):j+S
                        for ktt = 1:ifelse(pnrank==2,2,1)
                            if pnrank==2
                                if (S==0 && ktt==1) || (S==1 && ktt==2);itt=1;end
                                if (S==0 && ktt==2) || (S==1 && ktt==1);itt=0;end
                            end                    
                            if (-1)^ll1 != (-1)^ll2;continue;end
                            if itt != Int((1+(-1)^(ll1+S))/2);continue;end
                            if (-1)^(ll1+S+itt) != -1;continue;end
                            num=num+1
                            push!(pw_channels,[Tz,ll1,ll2,S,j])
                            dict_pwch[pnrank][[Tz,ll1,ll2,S,j]] = num
                            arr_pwch[pnrank][S+1][j+1][ll1-abs(j-S)+1][ll2-abs(j-S)+1] = num
                        end
                    end
                end
            end
        end
    end
    println(io,"# of two-body states $num")
    return pw_channels,dict_pwch,arr_pwch
end

"""
calculating coulomb contribution in HO base. Results are overwritten to `Vcoulomb`
"""
function calc_coulomb(chiEFTobj,Vcoulomb,pw_channels,num2bch;meshp=100)
    hw = chiEFTobj.hw
    Nnmax = chiEFTobj.Nnmax
    rmass = 0.5 * Mp
    brange = sqrt(hc2 / (rmass*hw))
    rmax = 20.0
    r,w = Gauss_Legendre(0.0,rmax,meshp)
    ra1 = zeros(Float64,meshp);rb1 = zeros(Float64,meshp)
    ra2 = zeros(Float64,meshp);rb2 = zeros(Float64,meshp)
    rnl1 = zeros(Float64,meshp);rnl2 = zeros(Float64,meshp)
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

"""
    Vrel(chiEFTobj::ChiralEFTobject,to) 
To define V in two-particle HO basis.
"""
function Vrel(chiEFTobj::ChiralEFTobject,to;calc_2n3n=false) 
    V12mom = ifelse(calc_2n3n,chiEFTobj.V12mom_2n3n,chiEFTobj.V12mom)
    pw_channels = chiEFTobj.pw_channels; Rnl = chiEFTobj.Rnl
    xr_fm = chiEFTobj.xr_fm; wr = chiEFTobj.wr; n_mesh = chiEFTobj.params.n_mesh
    Nnmax= chiEFTobj.params.Nnmax; num2bch = length(pw_channels)
    V12ab = [zeros(Float64,Nnmax+1,Nnmax+1) for i=1:num2bch]
    Vcoulomb = [zeros(Float64,Nnmax+1,Nnmax+1) for i=1:num2bch]
    if chiEFTobj.params.coulomb && !calc_2n3n 
        calc_coulomb(chiEFTobj.params,Vcoulomb,pw_channels,num2bch)
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
                Vab[n1+1,n2+1]= phase * vsum + ifelse(!calc_2n3n,t_vcoul,0.0)
            end
        end 
    end
    return V12ab
end

"""
    Calc_Deuteron(chiEFTobj::ChiralEFTobject,to;io=stdout)
Function to calculate deuteron binding energy.
"""
function Calc_Deuteron(chiEFTobj::ChiralEFTobject,to;io=stdout)
    V12mom = chiEFTobj.V12mom; pw_channels = chiEFTobj.pw_channels
    xr_fm = chiEFTobj.xr_fm; wr = chiEFTobj.wr; n_mesh = chiEFTobj.params.n_mesh
    ndim = 2*n_mesh
    H_d = zeros(Float64,ndim,ndim); V_d = zeros(Float64,ndim,ndim); T_d = zeros(Float64,ndim,ndim)
    ofst_i = ofst_j = 0
    Rmass = Mp * Mn / (Mp+Mn)
    for num in eachindex(V12mom)
        iz,l1,l2,isz,jz = pw_channels[num]
        if iz != 0 || isz !=1 || jz != 1; continue;end
        if l1 == l2 == 0
            ofst_i = ofst_j = 0
        elseif l1 == l2 == 2
            ofst_i = ofst_j = n_mesh
        elseif l1 == 0 && l2 == 2
            ofst_i = 0 ; ofst_j = n_mesh
        else
            continue
        end
        Vtmp = V12mom[num]
        for i = 1:n_mesh
            for j = 1:n_mesh
                idx_i = i + ofst_i
                idx_j = j + ofst_j
                v = Vtmp[i,j] * xr_fm[i]*xr_fm[j] * sqrt(wr[i]*wr[j])
                V_d[idx_i,idx_j] = V_d[idx_j,idx_i] = v
            end
            T_d[i,i] = T_d[i+n_mesh,i+n_mesh] = (xr_fm[i]*hc)^2 / (2*Rmass)
        end
    end
    H_d .= V_d
    H_d .+= T_d
    evals,evecs = eigen(H_d)
    E_d = minimum(evals)
    #println(io,"Deuteron energy:",@sprintf("%12.6f",E_d)," MeV")
    return E_d
end

"""
    hw_formula(A,fnum)

empirical formula for harmonis oscillator parameter hw by mass number A
fnum=2: for sd-shell, Ref. J. Blomqvist and A. Molinari, Nucl. Phys. A106, 545 (1968).
"""
function hw_formula(A::Int,fnum::Int)
    hw = 0.0
    if fnum == 2        
        hw = 45.0 * (A^(-1.0/3.0)) -25.0 * (A^(-2.0/3.0))
    else
        hw = 41.0 * (A^(-1.0/3.0))
    end
    return hw
end
