"""
Reference: 
"""
function calc_nuclearmatter_in_momspace(chiEFTobj::ChiralEFTobject,to,io)
    @assert chiEFTobj.params.n_mesh >= 100 "n_mesh >=100 is recomended for E/A calc."
    n_below_kF = get_n_below_kF(chiEFTobj)
    kF = chiEFTobj.params.kF
    rho = (2 *kF^3 )/ (3*pi^2)
    println(io,"kF (fm^-1) $kF rho (fm^-3) $rho n_below_kF $n_below_kF")
    TperA = 3*hc^2 *kF^2 /(5*2* (Mp+Mn)/2)
    EHF  = calc_EperA_HF(chiEFTobj,n_below_kF)
    EPT2 = calc_EperA_PT2(chiEFTobj,to)
    E = TperA + EHF + EPT2
    println(io,"kinetic ",@sprintf("%10.3f",TperA)," EHF/A ",@sprintf("%10.3f",EHF)," EPT2/A ",@sprintf("%15.8f",EPT2),
    "   E(1) ", @sprintf("%10.3f",E-EPT2)," E(1+2) ",@sprintf("%10.3f",E))
    return nothing
end

function calc_EperA_HF(chiEFTobj,n_below_kF;verbose=false)
    V12mom = chiEFTobj.V12mom
    V12mom_2n3n = chiEFTobj.V12mom_2n3n
    pw_channels = chiEFTobj.pw_channels
    xr_fm = chiEFTobj.xr_fm
    wr = chiEFTobj.wr
    nthre = nthreads()
    sumarrs = zeros(Float64,nthre)
    @threads for chidx in eachindex(V12mom)
        tV = V12mom[chidx]; tV_2n3n = V12mom_2n3n[chidx]
        Tz,l1,l2,S,J = pw_channels[chidx]
        if l1 != l2; continue;end
        Jfac = 2*J + 1
        LSTfac = 1 - (-1)^(l1+S+1)
        tsum = 0.0
        afac = ifelse(Tz==0,sqrt(2.0),LSTfac/sqrt(2.0))
        for i = 1:n_below_kF
	        k = xr_fm[i]; k2 = k^2
            w = wr[i]
            v2 = tV[i,i]; v2n3n = tV_2n3n[i,i]
            tsum  += k2 * (1-1.5*(k/kF)+ 0.5 * (k/kF)^3 ) * w * Jfac *  afac^2 * (v2 + v2n3n/2)
        end
        sumarrs[threadid()] += tsum
        if verbose && abs(tsum) >  1.0                                                                        
           println("to HF: L $l1 S $S J $J tsum ",@sprintf("%10.4e",tsum))
        end
    end
    EHF = sum(sumarrs)
end
function calc_EperA_PT2(chiEFTobj::ChiralEFTobject,to)
    kF = chiEFTobj.params.kF
    n_mesh = chiEFTobj.params.n_mesh
    V12mom = chiEFTobj.V12mom
    V12mom_2n3n = chiEFTobj.V12mom_2n3n
    pw_channels = chiEFTobj.pw_channels
    xr_fm = chiEFTobj.xr_fm
    wr = chiEFTobj.wr
    nKmesh = 50
    Ks,wKs = Gauss_Legendre(0,2*kF,nKmesh)
    EPT2 = 0.0
    for idx1 in eachindex(V12mom)
        Tz,tLp,tL,S,J = pw_channels[idx1]
        v1Mat = V12mom[idx1]; v1Mat_2n3n = V12mom_2n3n[idx1]
        LSfac1 = 1 - (-1)^(tL+S+1)   
        LSfac2 = 1 - (-1)^(tLp+S+1)   
        afac1 = ifelse(Tz==0,sqrt(2.0),LSfac1/sqrt(2.0))
        afac2 = ifelse(Tz==0,sqrt(2.0),LSfac2/sqrt(2.0))
        afac = afac1 * afac2
        coeff = (2*J+1) / (16*kF^3)
        for iK = 1:nKmesh
            K = Ks[iK]; K2 = K^2
            Kfac = K2 * wKs[iK] 
            Ksq = sqrt(kF^2 - K^2/4)
            for i = 1:n_mesh
                ki = xr_fm[i]; wi = wr[i]
                if ki > Ksq;continue;end
                for j = 1:n_mesh                    
                    kj = xr_fm[j]; wj = wr[j]
                    if kj <= Ksq;continue;end
                    v1 = afac * (v1Mat[i,j] + v1Mat_2n3n[i,j])
                    momfac = (ki^2) *(kj^2)* wi *wj
                    edeno = get_spe_denominator(Tz,ki,kj,K)
                    EPT2 += coeff * momfac/edeno * v1^2 * Kfac
                end 
            end
        end
    end
    return EPT2
end

function get_spe_denominator(Tz,ki,kj,K,Lzero=true)
    mN = (Mp+Mn)/2
    mass = mN       
    if Tz == -2; mass = Mp; end
    if Tz == 2; mass = Mn; end
    return hc^2 * ((ki+0.5*K)^2 + (-ki+0.5*K)^2 - (kj+0.5*K)^2 - (-kj+0.5*K)^2) / mass
end

function calc_EperA_PT2(chiEFTobj::ChiralEFTobject,n_below_kF::Int,to)
    barLmax = 6
    kF = chiEFTobj.params.kF
    mN = (Mp + Mn)/2
    n_mesh = chiEFTobj.params.n_mesh
    d6j_int = chiEFTobj.d6j_int
    V12mom = chiEFTobj.V12mom
    V12mom_2n3n = chiEFTobj.V12mom_2n3n
    pw_channels = chiEFTobj.pw_channels
    xr_fm = chiEFTobj.xr_fm
    wr = chiEFTobj.wr
    numts = 24
    ts,ws = Gauss_Legendre(-1,1,numts)
    PLs = [ zeros(Float64,numts) for L=0:barLmax] 
    for L = 0:barLmax
        for i in eachindex(ts)
            PLs[L+1][i] = Legendre(L,ts[i])
        end
    end
    nKmesh = 50
    Ks,wKs = Gauss_Legendre(0,2*kF,nKmesh)
    EPT2 = 0.0
    nthre = nthreads()
    sumarrs = zeros(Float64,nthre)
    @threads for idx1 in eachindex(V12mom)
        Tz,tLp,tL,S,tJ = pw_channels[idx1]
        v1Mat = V12mom[idx1]; v1Mat_2n3n = V12mom_2n3n[idx1]
        sqfac1 = hat(tL)*hat(tLp)
        LSfac1 = 1 - (-1)^(tL+S+1)   
        tid = threadid()
        if LSfac1 == 0;continue;end
        for idx2 in eachindex(V12mom)
            Tz2,Lp,L,S2,J = pw_channels[idx1]
            if S != S2; continue;end
            LSfac2 = 1 - (-1)^(L+S+1)
            if LSfac2 == 0;continue;end
            v2Mat = V12mom[idx2]; v2Mat_2n3n = V12mom_2n3n[idx2]
            sqfac2 = hat(L)*hat(Lp)
            Jfac = (2*J+1) *(2*tJ+1)
            for barL = 0:barLmax
                PL = PLs[barL+1]
                if !tri_check(L,tLp,barL);continue;end
                if !tri_check(Lp,tL,barL);continue;end
                if !tri_check(tJ,barL,J);continue;end
                if !tri_check(L,barL,tLp);continue;end
                if !tri_check(J,barL,tJ);continue;end
                if !tri_check(tL,barL,Lp);continue;end
                ifac = (-1)^(tL+Lp+barL) * (-1)^div(L-Lp+tL-tLp,2) #added
                wsym1 = clebschgordan(Float64,L,0,tLp,0,barL,0) *get_dict6jint_for9j(L,S,J,tJ,barL,tLp,d6j_int)
                wsym2 = clebschgordan(Float64,Lp,0,tL,0,barL,0) *get_dict6jint_for9j(J,S,Lp,tL,barL,tJ,d6j_int)
                coeff = ifac * sqfac1 * sqfac2 *Jfac * LSfac1 * LSfac2 * wsym1 *wsym2 * (4*pi) 
                if coeff == 0.0;continue;end
                for iK = 1:nKmesh
                    K = Ks[iK]; K2 = K^2
                    Kfac = K2 * wKs[iK] 
                    Ksq = sqrt(kF^2 - K^2/4)
                    for i = 1:n_mesh
                        ki = xr_fm[i]; wi = wr[i]
                        if ki > Ksq;continue;end
                        for j = 1:n_mesh                    
                            kj = xr_fm[j]; wj = wr[j]
                            if kj <= Ksq;continue;end
                            v1 = v1Mat[j,i] + v1Mat_2n3n[j,i]/2
                            v2 = v2Mat[i,j] + v2Mat_2n3n[i,j]/2
                            momfac = (ki^2) *(kj^2)* wi *wj
                            edeno = get_spe_denominator(Tz,ki,kj,K)
                            for idx_t in eachindex(ts) # theta_pp' integral
                                Pw = PL[idx_t] * ws[idx_t]
                                sumarrs[threadid()] += coeff * momfac/edeno * v1 * v2 * Pw *Kfac
                            end
                        end
                    end
                end
            end
        end
    end
    EPT2 = sum(sumarrs)
    rho = (2 *kF^3 )/ (3*pi^2)
    return EPT2 /(4 * pi^5) 
end

function Legendre(n,x)
    nmax = Int(ifelse(n%2==0,n//2,(n-1)//2))
    tsum = 0.0
    for k = 0:nmax
        tsum += (-1)^k * factorial(big(2*n-2*k)) / factorial(k) / factorial(n-k) / factorial(n-2*k) *x^(n-2*k)
    end
    ret = tsum / 2^n
    return ret
end

function get_n_below_kF(chiEFTobj;verbose=true)
    kF = chiEFTobj.params.kF  
    xr_fm = chiEFTobj.xr_fm
    n_below_kF = 0
    for i = 1:length(xr_fm)
        if xr_fm[i] > kF; break; end
        n_below_kF = i
    end    
    return n_below_kF 
end
