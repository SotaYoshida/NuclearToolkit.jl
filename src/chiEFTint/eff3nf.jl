function XF(Anum)
    return (1.0 + (0.35*pi/1.15)^2 * Anum^(-2/3) )^(-1/3)
end
function prep_Fis!(chiEFTobj,xr,F0s,F1s,F2s,F3s)
    if !chiEFTobj.calc_3N;return nothing;end
    kF = chiEFTobj.kF
    n_mesh = chiEFTobj.n_mesh
    mpi = sum(mpis)/3.0 / hc; mpi2 = mpi^2; mpi4=mpi2^2
    kF2 = kF^2; kF4 = kF2^2
    fac = 1.0 / (2*pi)^2
    @threads for i = 1:n_mesh
        k = xr[i] / hc; k2 = k^2; k4=k2^2
        F0s[i] = fac * ( kF
                         + (kF2+mpi2-k2)/(4.0*k) *log( ((k+kF)^2 +mpi2)/((k-kF)^2 +mpi2))
                         - mpi * (atan((k+kF)/mpi)-atan((k-kF)/mpi)))
        F1s[i] = fac * ( kF/(4.0*k2) *(3.0*k2-kF2-mpi2)
                         - mpi * (atan((k+kF)/mpi)-atan((k-kF)/mpi))
                         + 1.0/(16.0 * k^3) * ( mpi4 +2*mpi2*(3.0*k2+kF2)
                                                +(kF2-k2)*(kF2+3.0*k2))*log(((k+kF)^2 +mpi2)/((k-kF)^2 +mpi2)))
        F2s[i] = fac/k2 * (kF*(3.0*k2+kF2-9.0*mpi2)/6.0
                           + (kF4-k4-mpi4+6.0*k2*mpi2)/(8.0*k) * log(((k+kF)^2 +mpi2)/((k-kF)^2 +mpi2))
                           +mpi*(mpi2-k2)*(atan((k+kF)/mpi)-atan((k-kF)/mpi)))
        F3s[i] = fac/k2 * ( ( (kF2+k2+mpi2)^3 + 2.0*k2*mpi2*(mpi2+6.0*k2)
                              -2.0*k2*(kF4+3.0*k4)) /(32.0*k^3) * log(((k+kF)^2 +mpi2)/((k-kF)^2 +mpi2))
                            -k2 * mpi * (atan((k+kF)/mpi)-atan((k-kF)/mpi))
                            -kF/(8.0*k2) * (mpi4 + 4.0*k2*mpi2 + kF4 -5.0*k4 +4.0*k2*kF2/3.0 + 2.0*mpi2*kF2))
    end
    return nothing
end

"""
    prep_integrals_for2n3n(chiEFTobj)

preparing integrals and Wigner symbols for density-dependent 3NFs:
- `Fis::Fis_2n3n` vectors {F0s,F1s,F2s,F3s}, Eqs.(A13)-(A16) in [Kohno2013].
- `QWs::` second kind Legendre functions, Eqs.(B1)-(B5) in [Kohno2013].
- `wsyms::` Wingner symbols with specific `j` used for both 2n and 2n3n.

Reference:  
[Kohno2013] M.Kohno Phys. Rev. C 88, 064005(2013)
"""
function prep_integrals_for2n3n(chiEFTobj,xr,ts,ws)
    calc_3N = chiEFTobj.calc_3N
    n_mesh = ifelse(calc_3N,chiEFTobj.n_mesh,3)
    F0s = zeros(Float64,n_mesh); F1s = zeros(Float64,n_mesh); F2s = zeros(Float64,n_mesh); F3s = zeros(Float64,n_mesh)
    QWs = prep_QWs(chiEFTobj,xr,ts,ws)
    prep_Fis!(chiEFTobj,xr,F0s,F1s,F2s,F3s)
    wsyms = prep_wsyms()
    Fis = Fis_2n3n(F0s,F1s,F2s,F3s)
    return util_2n3n(Fis,QWs,wsyms)
end

function calc_vmom_3nf(chiEFTobj,it,to;pnm=false)
    if !chiEFTobj.params.calc_3N; return nothing; end
    if it > 1; for i in eachindex(chiEFTobj.pw_channels); chiEFTobj.V12mom_2n3n[i] .= 0.0; end;end
    dLECs = chiEFTobj.LECs.dLECs; xr = chiEFTobj.xr
    V12mom = chiEFTobj.V12mom; dict_pwch = chiEFTobj.dict_pwch
    util_2n3n = chiEFTobj.util_2n3n; lsjs = chiEFTobj.lsjs
    tmp_llsj = chiEFTobj.tllsj                    
    kF = chiEFTobj.params.kF
    n_mesh = chiEFTobj.params.n_mesh
    LamChi = 700.0/hc
    rho = ifelse(pnm,2.0,4.0) * kF^3 / (6.0 * pi^2)
    nreg = 3
    mpi = sum(mpis)/3.0/hc; mpi2 = mpi^2
    bbf =  2.0/3.0 / (4.0 * pi^2)

    r_c1 = dLECs["ct1_NNLO"] *1.e-3 *gA^2 * mpi2  / ((Fpi/hc)^4) *hc^2
    r_c3 = dLECs["ct3_NNLO"] *1.e-3 *gA^2 / (2.0*(Fpi/hc)^4) *hc^2
    r_c4 = dLECs["ct4_NNLO"] *1.e-3 *gA^2 / (2.0*(Fpi/hc)^4) *hc^2
    r_cD = dLECs["cD"] * hc * abs(gA) / (8.0*(Fpi/hc)^4 *LamChi)
    r_cE= -6.0/4.0 * dLECs["cE"] * rho / ( (Fpi/hc)^4 *LamChi) * hc
    # Contact term
    V_E(chiEFTobj.params,xr,r_cE*bbf,V12mom,dict_pwch,to)
    ## pion exchange: V_C(TPE), V_D(OPE)    
    tllsj = [copy(tmp_llsj) for i =1:nthreads()]
    for pnrank =1:3
        tdict = dict_pwch[pnrank]
        #MN = Ms[pnrank]#;dwn = 1.0/MN;sq_dwn=dwn^2      
        itt = 2 *(pnrank -2)
        @views tllsj[1:nthreads()][1] .= itt
        @inbounds for J=0:jmax
            lsj = lsjs[J+1]
            @inbounds for i= 1:n_mesh
                x = xr[i]
                @inbounds for j = 1:n_mesh
                    y = xr[j]
                    fac3n = bbf  *freg(x,y,nreg)
                    ## eff3nf 
                    single_3NF_pe(J,pnrank,x,y,mpi2,fac3n,rho,
                                  r_c1,r_c3,r_c4,r_cD,r_cE,
                                  util_2n3n,lsj,tllsj[threadid()],
                                  tdict,V12mom,i,j,to)
                end
            end
        end
    end
    return nothing
end

function V_E(chiEFTobj,xr,cE,V12mom,dict_pwch,to)
    nreg = 3
    pfunc = f1
    for pnrank = 1:3
        l=0;lp=0;S=0;J=0 ##1S0
        calc_Vmom!(chiEFTobj,pnrank,V12mom,dict_pwch[pnrank],xr,cE,cE,
                   l,lp,S,J,pfunc,nreg,to;is_3nf=true)
        l=0;lp=0;S=1;J=1 ##3S1
        if pnrank%2==1;continue;end
        calc_Vmom!(chiEFTobj,pnrank,V12mom,dict_pwch[pnrank],xr,cE,cE,
                   l,lp,S,J,pfunc,nreg,to;is_3nf=true)
    end
    return nothing
end

function S12(ell,ellp,J)
    r = 0.0
    if abs(ell-ellp)==2
        r = 6.0*sqrt(J*(J+1)) /(2*J+1)
    else
        if ell == ellp == J
            r = 2.0
        elseif ell==ellp== J+1
            r = -2.0*(J+2) / (2*J+1)
        elseif ell==ellp==J-1
            r = -2.0*(J-1) / (2*J+1)
        end
    end
    return r
end

""" 
    prep_QWs(chiEFTobj,xr,ts,ws)

returns struct `QLs`, second kind Legendre functions, Eqs.(B1)-(B5) in [Kohno2013].
Note that QLs.QLdict is also used in OPEP to avoid redundant calculations.
Reference:  [Kohno2013] M.Kohno Phys. Rev. C 88, 064005(2013)
"""
function prep_QWs(chiEFTobj,xr,ts,ws)
    n_mesh = chiEFTobj.n_mesh
    kF = chiEFTobj.kF
    dim = lmax + 2
    mpi = sum(mpis)/3.0 / hc; mpi2 = mpi^2
    QL0s  = [ zeros(Float64,n_mesh,n_mesh) for ell =0:dim]
    QWL2s = [ zeros(Float64,n_mesh,n_mesh) for ell =0:dim]
    QWL1s = [ zeros(Float64,n_mesh,n_mesh) for ell =0:dim]
    QWL0s = [ zeros(Float64,n_mesh,n_mesh) for ell =0:dim]
    QWlss = [ zeros(Float64,n_mesh,n_mesh) for ell = 0:dim]
    ndQW1s = [ [ zeros(Float64,n_mesh,n_mesh) for ellp = 0:dim] for ell=0:dim]
    ndQW2s = [ [ zeros(Float64,n_mesh,n_mesh) for ellp = 0:dim] for ell=0:dim]
    QLdict = [ Dict{Float64,Float64}() for ell=0:dim+1]
    for i = 1:n_mesh #for (i,x) in enumerate(xr)
        x = xr[i]
        x_fm = x/hc
        for (j,y) in enumerate(xr)
            y_fm = y/hc
            z = (x^2 + y^2 + mpi2*hc^2) / (2.0 * x*y)
            for ell = 0:lmax+1
                QL0s[ell+1][i,j] = QL(z,ell,ts,ws,QLdict)
                QWL0s[ell+1][i,j] = QWL(kF,x_fm,y_fm,ell,ell,0,ts,ws,mpi2,QLdict)   
            end
            for ell=0:dim
                QWL2s[ell+1][i,j] = QWL(kF,x_fm,y_fm,ell,ell,2,ts,ws,mpi2,QLdict)
                QWL1s[ell+1][i,j] = QWL(kF,x_fm,y_fm,ell,ell,1,ts,ws,mpi2,QLdict)
                QWlss[ell+1][i,j] = QWls(kF,x_fm,y_fm,ell,ell,0,ts,ws,mpi2,QLdict)
                for ellp = 0:dim
                    ndQW1s[ell+1][ellp+1][i,j] = QWL(kF,x_fm,y_fm,ell,ellp,1,ts,ws,mpi2,QLdict;is_xmul=false)
                    ndQW2s[ell+1][ellp+1][i,j] = QWL(kF,x_fm,y_fm,ell,ellp,2,ts,ws,mpi2,QLdict;is_xmul=false)
                end 
            end
        end
    end    
    return QLs(QL0s,QWL2s,QWL1s,QWL0s,QWlss,ndQW1s,ndQW2s,QLdict)
end

function single_3NF_pe(J,pnrank,x,y,mpi2,fac3n,rho,
                       c1,c3,c4,cD,cE,util_2n3n,
                       lsj,tllsj,tdict,V12mom,V_i,V_j,to)
    allsum = 0.0
    f_idx = 6; if J==0;f_idx = 3;end
    F0s = util_2n3n.Fis.F0s; F1s = util_2n3n.Fis.F1s
    F2s = util_2n3n.Fis.F2s; F3s = util_2n3n.Fis.F3s
    QL0s = util_2n3n.QWs.QL0s ;  QWL2s = util_2n3n.QWs.QWL2s
    QWL1s = util_2n3n.QWs.QWL1s; QWL0s = util_2n3n.QWs.QWL0s; QWlss = util_2n3n.QWs.QWlss
    ndQW1s = util_2n3n.QWs.ndQW1s; ndQW2s = util_2n3n.QWs.ndQW2s
    wsyms = util_2n3n.wsyms
    z  = (x^2 + y^2 + mpi2*hc^2) / (2.0 * x*y)
    x_fm = x/hc; x_fm2 = x_fm^2
    y_fm = y/hc; y_fm2 = y_fm^2
    XYT = 2.0*x_fm*y_fm

    QWJ0p1 = QWL0s[J+2][V_j,V_i]
    QWJ0m1 = 0.0;if J>0;QWJ0m1 = QWL0s[J][V_j,V_i];end
    QJ   = QL0s[J+1][V_i,V_j]
    QJp1 = QL0s[J+1+1][V_i,V_j]        
    QJ1 = z*QJ - delta(J,0)        
    QDL_J = -(J+1)/(z^2 -1.0) * (z*QJ - QJp1)        
    QWL2_J  = QWL2s[J+1][V_j,V_i]    
    F0X = F0s[V_i]; F1X = F1s[V_i]; F2X = F2s[V_i]; F3X = F3s[V_i]
    F0Y = F0s[V_j]; F1Y = F1s[V_j]; F2Y = F2s[V_j]; F3Y = F3s[V_j]
    
    @inbounds for idx = 1:f_idx
        @views tllsj[2:5] .= lsj[idx]
        V12idx = get(tdict,tllsj,-1)
        if V12idx == -1;continue;end
        allsum = 0.0
        tV12 = V12mom[V12idx]        
        ell,ellp,S,Jrel = lsj[idx] #Jrel=J
        if pnrank%2==1 && (ell+S) %2 ==1;continue;end
        QL0  = QL0s[ell+1][V_i,V_j]
        QL0p = QL0s[ellp+1][V_i,V_j]
        QLm1 = 0.0; if ell > 0; QLm1= QL0s[ell][V_i,V_j]; end
        QL1 = z*QL0 - delta(ell,0) # Ql^(1)=zQl -delta(l,0)
        QL1p = z*QL0p - delta(ellp,0)
        QLp1  = QL0s[ell+2][V_i,V_j]        
        QLp1p = QL0s[ellp+2][V_i,V_j]
        QLp2 = QL0s[ell+2+1][V_i,V_j]

        QDL_ell = -(ell+1)/(z^2 -1.0) *(z*QL0 - QLp1) 
        QDL_ellp = -(ellp+1)/(z^2 -1.0) *(z*QL0p - QLp1p)
        QWL2_l  = QWL2s[ell+1][V_i,V_j]
        QWL2_lp = QWL2s[ellp+1][V_j,V_i]
        
        QW2_lm1 = 0.0
        if ell>=1; QW2_lm1 = QWL2s[ell][V_i,V_j];end
        QW2_lp1 = QWL2s[ell+2][V_i,V_j]
        QWL1_xyl = QWL1s[ell+1][V_i,V_j]
        QWL1_yxl = QWL1s[ell+1][V_j,V_i]
        QWL0xy = QWL0s[ell+1][V_i,V_j]
        QWL0yx = QWL0s[ell+1][V_j,V_i]

        QWL0_lp1 = QWL0s[ell+2][V_i,V_j]
        QWL0_lm1 = 0.0
        if ell-1 >= 0;QWL0_lm1 = QWL0s[ell][V_i,V_j]; end
        QWL0_lp2yx = 0.0
        if ell+3 <= lmax+1; QWL0_lp2yx = QWL0s[ell+3][V_i,V_j];end
        QWL0_lm2 = 0.0;if ell-2 >= 0; QWL0_lm2 = QWL0s[ell-1][V_i,V_j];end
        Wls0 = QWlss[ell+1][V_i,V_j]
        Wls0m1 = 0.0; if ell>=1;Wls0m1=QWlss[ell][V_i,V_j];end
        Wls0p1 = 0.0; if ell+2<=lmax+1;Wls0p1=QWlss[ell+2][V_i,V_j];end

        QX1 = QX(ell,ndQW1s,V_i,V_j) + QX(ell,ndQW1s,V_j,V_i)
        QX1p= QX(ellp,ndQW1s,V_i,V_j) + QX(ellp,ndQW1s,V_j,V_i)
        QX1J = QX(J,ndQW1s,V_i,V_j) + QX(J,ndQW1s,V_j,V_i)
        QX1_p1 = QX(ell+1,ndQW1s,V_i,V_j) + QX(ell+1,ndQW1s,V_j,V_i) 
        QX1_m1 = 0.0
        if ell >=2; QX1_m1= QX(ell-1,ndQW1s,V_i,V_j) + QX(ell-1,ndQW1s,V_j,V_i);end

        QDL_lm1 = -(ell)/(z^2 -1.0) *(z*QLm1 - QL0) 
        QLp2 = QL0s[ell+2+1][V_i,V_j]
        QDL_lp1 =  -(ell+2)/(z^2 -1.0) *(z*QLp1 - QLp2)  
        
        d3jinv = 1.0/l2l[ell+1]
        d3jinvnd = 1.0/l2lnd[ellp+1][ell+1]

        ttis = ifelse(S%2==0,-3.0,1.0)
        ttit = ifelse((S+ell)%2==1,-3.0,1.0)
        tS12 = S12(ell,ellp,J)

        ## Central term
        if ell == ellp
            ## V_D: central Eq.(B18)
            tmp = Vc_cD(ell,XYT,rho,mpi2,ttit,ttis,F0X,F0Y,QL0)
            allsum += tmp * fac3n * cD
            ## V_c1: central Eq.(B6)
            tmp = Vc_c1(ell,x_fm2,y_fm2,XYT,rho,mpi2,ttit,ttis,
                        F0X,F0Y,F1X,F1Y,QL0,QWL0_lp1,QWL0_lm1,
                        QL1,QWL1_xyl,QWL1_yxl,QWL2_l,QDL_ell)
            allsum += tmp * fac3n * c1
            ## V_c3: central Eq.(B11)
            tmp = Vc_c3(ell,x_fm,y_fm,x_fm2,y_fm2,XYT,rho,mpi2,ttit,ttis,
                        F0X,F0Y,F1X,F1Y,F2X,F2Y,F3X,F3Y,
                        QL0,QWL0yx,QWL0_lp1,QWL0_lm1,QWL0_lp2yx,QWL0_lm2,
                        QDL_ell)
            allsum +=  tmp * fac3n * c3
            ## V_c4: central Eq.(B15)
            tmp = Vc_c4(ell,x_fm2,y_fm2,XYT,rho,mpi2,ttit,ttis,
                        F0X,F0Y,F1X,F1Y,F2X,F2Y,F3X,F3Y,
                        QL0,QWL0yx,QWL0_lp1,QWL0_lm1,QWL0_lp2yx,QWL0_lm2) 
            allsum += tmp * fac3n * c4
        end
        ## Tensor term
        if (ellp == ell+2 && J==ell+1) || (ellp == ell-2 && J==ell-1)
            ## V_c1: tensor non-diagonal Eq.(B7)
            tmp = Vt_c1_nd(x_fm2,y_fm2,XYT,rho,ttit,ttis,tS12,F0X,F0Y,F1X,F1Y,
                           QL0,QL0p,QJ,QDL_ell,QDL_ellp,QDL_J)
            allsum += tmp * fac3n * c1
            ## V_c3: tensor non-diagonal Eq.(B11)
            tmp =Vt_c3_nd(x_fm2,y_fm2,XYT,rho,mpi2,ttit,ttis,tS12,
                          F0X,F0Y,F1X,F1Y,F2X,F2Y,F3X,F3Y,
                          QL0,QL0p,QJ,QL1,QL1p,QJ1,QDL_ell,QDL_ellp,QDL_J)
            allsum += tmp * fac3n * c3
            ## V_c4: tensor non-diagonal Eq.(B16+erratum)
            tmp = Vt_c4_nd(V_i,V_j,ell,ellp,J,x_fm2,y_fm2,XYT,rho,mpi2,
                           F0X,F0Y,F1X,F1Y,F2X,F2Y,F3X,F3Y,QL0,QL0p,QJ,
                           QWL0yx,QWL0_lp1,QWL0_lm1,QWL0_lp2yx,QWL0_lm2,
                           QWJ0p1,QWJ0m1,QWL2_l,QWL2_lp,QWL2_J,
                           QX1,QX1p,QX1J,ndQW1s,ndQW2s,wsyms,d3jinvnd)
            allsum += tmp * tS12 * ttit * fac3n * c4
            ## V_D: tensor non-diagonal,  Eq.(B19)
            tmp =Vt_cD_nd(ell,ellp,x_fm2,y_fm2,XYT,rho,F0X,F0Y,F1X,F1Y,F3X,F3Y,
                          QL0,QL0p,QJ,tS12,ttit) 
            allsum += tmp * fac3n * cD
        end
        if ellp == ell >=1 && abs(ell-J) <= 1 && S==1
            ## V_c1: tensor diagonal Eq.(B8)
            tmp = Vc_c1(ell,x_fm2,y_fm2,XYT,rho,ttit,ttis,tS12,F0X,F0Y,F1X,F1Y,
                        QL0,QL0p,QLm1,QLp1,QDL_ell,QDL_lm1,QDL_lp1)
            allsum += tmp * fac3n * c1
            ## V_c3: tensor diagonal Eq.(B12)
            tmp= Vt_c3_d(V_i,V_j,ell,ellp,J,x_fm2,y_fm2,XYT,z,rho,mpi2,ttit,tS12,
                         F0X,F0Y,F1X,F1Y,F2X,F2Y,F3X,F3Y,QL0,QL0p,QJ,QLm1,QLp1,
                         QL1,QWL0yx,QWL0_lp1,QWL0_lm1,QWL0_lp2yx,QWL0_lm2,
                         QWJ0p1,QWJ0m1,QWL2_l,QWL2_lp,QWL2_J,QW2_lm1,QW2_lp1,
                         QX1,QX1p,QX1J,QX1_m1,QX1_p1,QDL_ell,QDL_lm1,QDL_lp1,
                         ndQW1s,ndQW2s,wsyms,d3jinv)
            allsum += tmp * fac3n * c3
            ## V_c4: tensor diagonal Eq.(B17+Erratum) 
            tmp = Vt_c4_d(V_i,V_j,ell,ellp,J,x_fm2,y_fm2,XYT,rho,mpi2,
                          F0X,F0Y,F1X,F1Y,F2X,F2Y,F3X,F3Y,
                          QL0,QL0p,QJ,QLm1,QLp1,
                          QWL0yx,QWL0_lp1,QWL0_lm1,QWL0_lp2yx,QWL0_lm2,
                          QWJ0p1,QWJ0m1,QWL2_l,QWL2_lp,QWL2_J,QW2_lm1,QW2_lp1,
                          QX1,QX1p,QX1J,QX1_m1,QX1_p1,ndQW1s,ndQW2s,wsyms,d3jinv)
            allsum += tmp * tS12 * ttit * fac3n * c4
            ## V_D: tensor diagonal Eq.(B20)
            tmp = Vt_cD_d(ell,x_fm2,y_fm2,XYT,rho,ttit,ttis,tS12,
                          F0X,F0Y,F1X,F1Y,F3X,F3Y,QL0,QL0p,QLm1,QLp1)
            allsum += tmp * fac3n * cD
        end
        ## LS term
        if S==1 && ell == ellp >0
            ## V_c1: spin-orbit  Eq.(B9)
            facL = 3.0 * (ell*(ell+1)+2-J*(J+1))/(2*ell+1)
            QT = QWL0_lp1 - QWL0_lm1 +Wls0
            allsum +=  QT * fac3n *facL * c1
            ## V_c3 + V_c4: spin-orbit Eq.(B14)+erratum
            facL = 3.0 * (ell*(ell+1)+2-J*(J+1))/(2*ell+1)
            QLS  = (mpi2+0.5*(x_fm2+y_fm2)) * (QWL0_lm1 - QWL0_lp1 -Wls0)
            QLS += 0.5*XYT * ((ell-1)/(2*ell-1) *Wls0m1 + (ell+2)/(2*ell+3)*Wls0p1)
            if ell == 1; QLS += - XYT/4.0 * (F0X+F0Y-F1X-F1Y); end 
            allsum +=  QLS * fac3n *facL * (c3 + c4*ttit/3.0)
        end
        tV12[V_i,V_j] += allsum
    end
    return nothing
end

## V_c1: central part  Eq.(B6)
function Vc_c1(ell,x_fm2,y_fm2,XYT,rho,mpi2,ttit,ttis,
               F0X,F0Y,F1X,F1Y,QL0,QWL0_lp1,QWL0_lm1,
               QL1,QWL1_xyl,QWL1_yxl,QWL2_l,QDL_ell)
    QA = (QL0 + mpi2*QDL_ell/XYT) / XYT            
    QB = QL0/XYT * (x_fm2 * (F0X-F1X) + y_fm2 * (F0Y-F1Y)) 
    QB += -0.5* QL1 * (F0X+F0Y-F1X-F1Y)
    QC = QWL2_l - QWL1_xyl -QWL1_yxl
    QC += (ell+1)/(2*ell+1) *QWL0_lp1 + ifelse(ell>=1,ell/(2*ell+1) * QWL0_lm1,0.0)
    if ell == 0;QC += 0.5 * (F0X+F0Y); end
    return (rho*QA-2.0*QB)*ttis*ttit/3.0 -6.0*QC    
end

## V_c1: tensor term non-diagonal Eq.(B7)
function Vt_c1_nd(x_fm2,y_fm2,XYT,rho,ttit,ttis,tS12,F0X,F0Y,F1X,F1Y,
                  QL0,QL0p,QJ,QDL_ell,QDL_ellp,QDL_J)
    Q1TA  = -1.0/XYT * ( (y_fm2 * QDL_ell +x_fm2*QDL_ellp)/XYT -QDL_J)
    Q1TB  = x_fm2/XYT * QL0p * (F0Y-F1Y) 
    Q1TB += y_fm2/XYT * QL0 * (F0X-F1X)
    Q1TB += -0.5 * QJ * (F0X+F0Y-F1X-F1Y)
    return (Q1TA*rho/3.0 -2.0/3.0 *Q1TB) * tS12 * ttit
end

## V_c1: tensor term diagonal Eq.(B8)
function Vc_c1(ell,x_fm2,y_fm2,XYT,rho,ttit,ttis,tS12,
               F0X,F0Y,F1X,F1Y,QL0,QL0p,QLm1,QLp1,QDL_ell,QDL_lm1,QDL_lp1)
    Q1TA = -1.0/XYT * ( (x_fm2 +y_fm2) * QDL_ell /XYT 
                        -0.5 * ((2*ell+3)/(2*ell+1)*QDL_lm1 + (2*ell-1)/(2*ell+1)*QDL_lp1))
    Q1TB1 = (x_fm2*QL0p *(F0Y-F1Y)
             +y_fm2*QL0  *(F0X-F1X))/XYT
    Q1TB2 = -0.25 *(F0X+F0Y-F1X-F1Y)
    Q1TB2 *= (2*ell+3)/(2*ell+1)*QLm1 + (2*ell-1)/(2*ell+1)*QLp1
    Q1TB = Q1TB1 +Q1TB2
    return (rho * Q1TA/3.0 -2.0/3.0 *Q1TB ) * tS12 * ttit
end

## V_c3: central part  Eq.(B11)
function Vc_c3(ell,x_fm,y_fm,x_fm2,y_fm2,XYT,rho,mpi2,ttit,ttis,
               F0X,F0Y,F1X,F1Y,F2X,F2Y,F3X,F3Y,
               QL0,QWL0yx,QWL0_lp1,QWL0_lm1,QWL0_lp2yx,QWL0_lm2,
               QDL_ell)
    C3R  = delta(ell,0) - 2.0*mpi2/XYT *QL0 - (mpi2/XYT)^2 * QDL_ell
    C3B  = 1.0/3.0 *(-mpi2/XYT *QL0) * (x_fm2*F2X+y_fm2 *F2Y)
    C3B += 0.25 *(y_fm2-x_fm2+mpi2)^2 /XYT *QL0 * (F0X-2.0*F1X)
    C3B += 0.25 *(x_fm2-y_fm2+mpi2)^2 /XYT *QL0 * (F0Y-2.0*F1Y)
    C3B += ((y_fm/XYT)^2 * (x_fm2+y_fm2+mpi2)^2 -y_fm2 -2.0/3.0 *mpi2) * QL0/XYT *x_fm2 *F3X
    C3B += ((x_fm/XYT)^2 * (x_fm2+y_fm2+mpi2)^2 -x_fm2 -2.0/3.0 *mpi2) * QL0/XYT *y_fm2 *F3Y
    C3G  = 1.0/(2.0*XYT) * (x_fm2+y_fm2+2.0*mpi2)^2 * QWL0yx 
    C3G += -(x_fm2+y_fm2+2.0*mpi2) * ( (ell+1)/(2*ell+1) *QWL0_lp1)            
    C3G += 0.5*XYT/(2*ell+1) * ( (ell^2 /(2*ell-1) +(ell+1)^2 /(2*ell+3)) * QWL0yx
                                 + (ell+1)*(ell+2)/(2*ell+3) * QWL0_lp2yx )
    if ell==0
        C3B += 1.0/3.0 * (x_fm2*F2X+y_fm2 *F2Y)
        C3B += -0.25 * (y_fm2-3.0*x_fm2+mpi2)*(F0X-2.0*F1X)
        C3B += -0.25 * (x_fm2-3.0*y_fm2+mpi2)*(F0Y-2.0*F1Y)
        C3B += -1.0/XYT * (y_fm2/XYT *(y_fm2+x_fm2+mpi2)-2.0/3.0 *XYT) *x_fm2 *F3X
        C3B += -1.0/XYT * (x_fm2/XYT *(y_fm2+x_fm2+mpi2)-2.0/3.0 *XYT) *y_fm2 *F3Y
        C3G += 0.125*rho 
        C3G += - (3.0*mpi2/4.0 +0.5*y_fm2+0.25*x_fm2)*F0Y
        C3G += - (3.0*mpi2/4.0 +0.5*x_fm2+0.25*y_fm2)*F0X
        C3G += 0.25 * (y_fm2 *F2Y + x_fm2 * F2X)
    end 
    if ell == 1
        C3B += -XYT/12.0 * ((F0X-2.0*F1X) + (F0Y-2.0*F1Y))
        C3B += -y_fm2/(3.0*XYT) * x_fm2 *F3X
        C3B += -x_fm2/(3.0*XYT) * y_fm2 *F3Y
        C3G += XYT/6.0 * (F0Y+F0X-0.5*(F1Y+F1X))
    end     
    if ell >= 1; C3G += - (x_fm2+y_fm2+2.0*mpi2) * ell/(2*ell+1) * QWL0_lm1; end 
    if ell >= 2; C3G += 0.5 *XYT *ell*(ell-1)/((2*ell+1)*(2*ell-1)) * QWL0_lm2;end
    return (rho/3.0 * C3R -2.0/3.0 *C3B)*ttis*ttit -6.0*C3G
end

## V_c3: tensor part, non-diagonal, Eq.(B11)
function Vt_c3_nd(x_fm2,y_fm2,XYT,rho,mpi2,ttit,ttis,tS12,
                  F0X,F0Y,F1X,F1Y,F2X,F2Y,F3X,F3Y,
                  QL0,QL0p,QJ,QL1,QL1p,QJ1,QDL_ell,QDL_ellp,QDL_J)
    Q3TA  = y_fm2/XYT *(QL0  + mpi2/XYT * QDL_ell) 
    Q3TA += x_fm2/XYT *(QL0p + mpi2/XYT * QDL_ellp) 
    Q3TA -= QJ + mpi2/XYT * QDL_J 
    Q3TB  = (F0X-2.0*F1X)*(-x_fm2*(0.5*QL1p-x_fm2/XYT *QL0p)
                           +0.5*XYT *(0.5*QJ1-x_fm2/XYT*QJ))
    Q3TB += (F0Y-2.0*F1Y)*(-y_fm2*(0.5*QL1-y_fm2/XYT *QL0)
                           +0.5*XYT *(0.5*QJ1-y_fm2/XYT*QJ))
    Q3TB += 1.0/3.0 * (x_fm2*F2X
                       +y_fm2*F2Y) * (y_fm2/XYT *QL0 + x_fm2/XYT *QL0p -QJ)
    Q3TB += 1.0/3.0 * (2*x_fm2*F3X-y_fm2*F3Y)*x_fm2/XYT *QL0p -0.5 *x_fm2*F3X * QL1p
    Q3TB += 1.0/3.0 * (2*y_fm2*F3Y-x_fm2*F3X)*y_fm2/XYT *QL0  -0.5 *y_fm2*F3Y * QL1
    Q3TB += -1.0/6.0 * (x_fm2*F3X+y_fm2*F3Y) *QJ
    Q3TB += XYT/4.0 * (F3Y+F3X)*QJ1
    QT = (Q3TA*rho/3.0 -2.0/3.0 *Q3TB)
    return QT * tS12 * ttit 
end

## V_c3: tensor part, diagonal, Eq.(B12)
function Vt_c3_d(V_i,V_j,ell,ellp,J,x_fm2,y_fm2,XYT,z,rho,mpi2,ttit,tS12,
                 F0X,F0Y,F1X,F1Y,F2X,F2Y,F3X,F3Y,QL0,QL0p,QJ,QLm1,QLp1,
                 QL1,QWL0yx,QWL0_lp1,QWL0_lm1,QWL0_lp2yx,QWL0_lm2,
                 QWJ0p1,QWJ0m1,QWL2_l,QWL2_lp,QWL2_J,QW2_lm1,QW2_lp1,
                 QX1,QX1p,QX1J,QX1_m1,QX1_p1,QDL_ell,QDL_lm1,QDL_lp1,
                 ndQW1s,ndQW2s,wsyms,d3jinv)
    
    Q3TA  = 1.0/XYT * (x_fm2 +y_fm2) * (QL0 + mpi2 * QDL_ell /XYT )
    Q3TA += -0.5 * ((2*ell+3)/(2*ell+1)*(QLm1 + mpi2*QDL_lm1/XYT)
                    + (2*ell-1)/(2*ell+1)*(QLp1 + mpi2*QDL_lp1/XYT))
    Q3TB  = (F0X-2.0*F1X)*(-x_fm2*(0.5*QL1-x_fm2/XYT *QL0)
                           +0.25*XYT *((2*ell+3)/(2*ell+1)*(0.5*(z*QLm1-delta(ell,1))
                                                             -x_fm2/XYT *QLm1)
                                       +(2*ell-1)/(2*ell+1)*(0.5*(z*QLp1)-
                                                             x_fm2/XYT *QLp1)))
    Q3TB += (F0Y-2.0*F1Y)*(-y_fm2*(0.5*QL1-y_fm2/XYT *QL0)
                           +0.25*XYT *( (2*ell+3)/(2*ell+1)*(
                               0.5*(z*QLm1-delta(ell,1))-y_fm2/XYT *QLm1)
                                        +(2*ell-1)/(2*ell+1)*(0.5*(z*QLp1)
                                                              -x_fm2/XYT *QLp1)))
    Q3TB += 1.0/3.0 * (x_fm2*F2X+y_fm2*F2Y) * (
        x_fm2/XYT * QL0p + y_fm2/XYT * QL0 
        - 0.5/(2*ell+1)* ((2*ell+3)*QLm1 + (2*ell-1)*QLp1))
    Q3TB += (2.0*x_fm2*F3X-y_fm2*F3Y)/3.0 *x_fm2/XYT *QL0 
    Q3TB += (2.0*y_fm2*F3Y-x_fm2*F3X)/3.0 *y_fm2/XYT *QL0 
    Q3TB += -0.5 * (x_fm2*F3X+y_fm2*F3Y) * QL1 
    Q3TB += -(x_fm2*F3X+y_fm2*F3Y)/12.0 * ((2*ell+3)*QLm1 + (2*ell-1)*QLp1)/(2*ell+1)
    Q3TB += 0.125*XYT *(F3X+F3Y) *((2*ell-1)*(z*QLp1)) /(2*ell+1) 
    if ell >= 1
        Q3TB += 0.125*XYT *(F3X+F3Y) *((2*ell+3)*(z*QLm1-delta(ell,1))) /(2*ell+1)
    end                     
    return  (rho * Q3TA/3.0 -2.0/3.0 *Q3TB ) * tS12 * ttit
end

## V_c4: central part  Eq.(B15)
function Vc_c4(ell,x_fm2,y_fm2,XYT,rho,mpi2,ttit,ttis,
               F0X,F0Y,F1X,F1Y,F2X,F2Y,F3X,F3Y,
               QL0,QWL0yx,QWL0_lp1,QWL0_lm1,QWL0_lp2yx,QWL0_lm2)
    QA = (0.5*rho -mpi2 *(F0X+F0Y)
          -(x_fm2*F2X+y_fm2*F2Y)/3.0) * (delta(ell,0)-mpi2/XYT *QL0)
    QA += (F0Y-2.0*F1Y)*(0.25*(x_fm2-3.0*y_fm2+mpi2)*delta(ell,0) 
                         + XYT/12.0 *delta(ell,1)
                         -0.25 * (x_fm2-y_fm2+mpi2)^2 *QL0/XYT )
    QA += (F0X-2.0*F1X)*(0.25*(y_fm2-3.0*x_fm2+mpi2)*delta(ell,0) 
                         +XYT/12.0 *delta(ell,1)
                         -0.25 * (y_fm2-x_fm2+mpi2)^2 *QL0/XYT )
    QA += F3Y * ( delta(ell,0) * (0.25 * (x_fm2+y_fm2+mpi2)-2.0*y_fm2/3.0) 
                  +XYT/12.0 * delta(ell,1)
                  -QL0/XYT * (0.25*(x_fm2+y_fm2+mpi2)^2 -x_fm2*y_fm2 - 2.0*mpi2*y_fm2/3.0))
    QA += F3X * ( delta(ell,0) * (0.25 * (x_fm2+y_fm2+mpi2)-2.0*x_fm2/3.0) 
                  +XYT/12.0 * delta(ell,1)
                  -QL0/XYT * (0.25*(x_fm2+y_fm2+mpi2)^2 -x_fm2*y_fm2 - 2.0*mpi2*x_fm2/3.0))
    QB  = delta(ell,0) * (0.125*rho  
                          +0.25 *(2.0*y_fm2+x_fm2-mpi2) *F0Y
                          +0.25 *(2.0*x_fm2+y_fm2-mpi2) *F0X
                          -0.25*(x_fm2*F2X+y_fm2*F2Y))
    QB += delta(ell,1) * 0.5 * XYT *((F1X+F1Y)/6.0 -(F0X+F0Y)/3.0)
    QB += -0.5 * (x_fm2+y_fm2)*(x_fm2+y_fm2+4.0*mpi2)*QWL0yx /XYT
    QB += (x_fm2+y_fm2+2.0*mpi2)/(2*ell+1) * ((ell+1)*QWL0_lp1 + ifelse(ell>0,ell*QWL0_lm1,0.0))
    QB += -0.5*XYT/(2*ell+1)*( (ell+2)*(ell+1)/(2*ell+3) *QWL0_lp2yx
                               + ((ell+1)^2 /(2*ell+3) + ell^2 / (2*ell-1))*QWL0yx 
                               + ifelse(ell>=2,(ell-1)*ell/(2*ell-1) *QWL0_lm2,0.0))                      
    QT =  2.0/3.0 *(QA+QB) *ttit*ttis
    return QT
end

## V_c4: tensor non-diagonal part  Eq.(B16+erratum)
function  Vt_c4_nd(V_i,V_j,ell,ellp,J,x_fm2,y_fm2,XYT,rho,mpi2,
                   F0X,F0Y,F1X,F1Y,F2X,F2Y,F3X,F3Y,QL0,QL0p,QJ,
                   QWL0yx,QWL0_lp1,QWL0_lm1,QWL0_lp2yx,QWL0_lm2,
                   QWJ0p1,QWJ0m1,QWL2_l,QWL2_lp,QWL2_J,
                   QX1,QX1p,QX1J,ndQW1s,ndQW2s,wsyms,d3jinvnd)
    cg1s = wsyms.cg1s;   cg2s = wsyms.cg2s
    d6_121 = wsyms.d6_121; d6_21  = wsyms.d6_21
    d6_222 = wsyms.d6_222; d9_12  = wsyms.d9_12

    Q4TA  = (0.5*rho-mpi2*(F0X+F0Y)-x_fm2*F2X/3.0-y_fm2*F2Y/3.0  )
    Q4TA *= (y_fm2*QL0+x_fm2*QL0p-XYT*QJ)/XYT 
    Q4TA += 0.5*XYT/(2.0*J+1.0) * (QWJ0p1 -QWJ0m1) 
    Q4TA += (F0Y-2.0*F1Y) * ( 
        y_fm2*(x_fm2-y_fm2+mpi2)*QL0/(2.0*XYT) 
        -0.25*(x_fm2-y_fm2+mpi2)*QJ 
        -0.5 *y_fm2 *delta(ell,0))
    Q4TA += (F0X-2.0*F1X) * ( 
        x_fm2*(y_fm2-x_fm2+mpi2)*QL0p/(2.0*XYT) 
        -0.25*(y_fm2-x_fm2+mpi2)*QJ 
        -0.5 * x_fm2 *delta(ellp,0))
    Q4TA += F3Y/3.0 * ( 
        y_fm2*(-y_fm2+3.0*x_fm2+3.0*mpi2)*QL0 /(2.0*XYT)
        -0.25 *(y_fm2+3.0*x_fm2+3.0*mpi2)*QJ +XYT/4.0 *QL0p
        -1.5*y_fm2 *delta(ell,0))
    Q4TA += F3X/3.0 * ( 
        x_fm2*(-x_fm2+3.0*y_fm2+3.0*mpi2)*QL0p /(2.0*XYT)
        -0.25 *(x_fm2+3.0*y_fm2+3.0*mpi2)*QJ +XYT/4.0 *QL0
        -1.5*x_fm2 *delta(ellp,0))
    ## term 5a           
    Q5a = 1.0/3.0 *(x_fm2 * (QX1 -delta(ell,0)*0.5*(F0X+F0Y))
                    +y_fm2 * (QX1p -delta(ellp,0)*0.5*(F0X+F0Y))- XYT* QX1J )
    ## term 5b 
    tcg = cg1s[ellp+1,J+1] * cg1s[ell+1,J+1]
    d6j = d6_121[ellp+1,J+1,ell+1]
    Q5b = sqrt(5.0/6.0) * (-1)^J *(x_fm2+y_fm2) * d3jinvnd * tcg * d6j *(
        ndQW1s[ellp+1][J+1][V_j,V_i] + ndQW1s[ell+1][J+1][V_i,V_j])
    ## term 5c&5d
    sum = 0.0
    for k=ell-1:2:ell+1
        if k <0;continue;end
        for jp = ellp-1:2:ellp+1
            if jp <0;continue;end
            for j = k-1:2:k+1
                if j <0;continue;end
                tcg = (-1)^(jp+1-j) * cg1s[jp+1,j+1]*cg1s[j+1,k+1] 
                tcg *= cg1s[ellp+1,jp+1] * cg1s[ell+1,k+1]
                if tcg == 0.0;continue;end
                d6j = d6_121[jp+1,j+1,k+1] * d6_21[jp+1,k+1,ell+1,ellp+1]
                if d6j == 0.0;continue;end
                sum += XYT * ndQW1s[jp+1][j+1][V_j,V_i] * hat(jp)*hat(j) * tcg *d6j
            end 
        end 
    end 
    Q5c = (-1)^J /sqrt(18.0) * d3jinvnd * sum
    sum = 0.0
    for k=ellp-1:2:ellp+1
        if k <0;continue;end
        for jp = ell-1:2:ell+1
            if jp <0;continue;end
            for j = k-1:2:k+1
                if j <0;continue;end
                tcg = (-1)^(jp+1-j) * cg1s[jp+1,j+1]*cg1s[j+1,k+1] 
                tcg *= cg1s[ell+1,jp+1] * cg1s[ellp+1,k+1]
                if tcg == 0.0;continue;end
                d6j = d6_121[jp+1,j+1,k+1] * d6_21[jp+1,k+1,ellp+1,ell+1]
                if d6j == 0.0;continue;end
                sum += XYT * ndQW1s[jp+1][j+1][V_i,V_j] * hat(jp)*hat(j) * tcg *d6j
            end 
        end 
    end 
    Q5d = (-1)^J /sqrt(18.0) * d3jinvnd * sum
    ## term 5e&5f
    sum = 0.0 
    for j = ell-1:2:ell+1
        if j <0;continue;end
        for jp = j-1:2:j+1
            if jp <0;continue;end
            tcg = (-1)^(jp+1-j) * cg1s[jp+1,j+1] * cg1s[ell+1,j+1] * cg2s[jp+1,ellp+1]
            d6j = d6_121[jp+1,j+1,ell+1] * d6_222[jp+1,ellp+1,ell+1]
            sum += y_fm2 * ndQW1s[jp+1][j+1][V_j,V_i] * (2*jp+1) * tcg * d6j                    
        end 
    end 
    Q5ef = sqrt(35.0/3.0) * (-1)^J * d3jinvnd * sum
    sum = 0.0
    for j = ellp-1:2:ellp+1
        if j <0;continue;end
        for jp = j-1:2:j+1
            if jp <0;continue;end
            tcg = (-1)^(jp+1-j) * cg1s[jp+1,j+1] * cg1s[ellp+1,j+1] * cg2s[jp+1,ell+1]
            d6j = d6_121[jp+1,j+1,ellp+1] * d6_222[jp+1,ell+1,ellp+1]
            sum += x_fm2 * ndQW1s[jp+1][j+1][V_i,V_j] * (2*jp+1) * tcg * d6j                    
        end 
    end 
    Q5ef += sqrt(35.0/3.0) * (-1)^J * d3jinvnd * sum
    ## term 5g & 5h 
    sum = 0.0
    for j = ellp-1:2:ellp+1
        if j <0;continue;end
        for k = j-1:2:j+1
            if k <0;continue;end
            tcg = (-1)^(ellp+1-j) * cg1s[ellp+1,j+1] * cg1s[k+1,j+1] * cg2s[k+1,ell+1]
            d6j = d6_121[ellp+1,j+1,k+1] * d6_222[k+1,ell+1,ellp+1]
            sum += x_fm2 * ndQW1s[ellp+1][j+1][V_j,V_i] * (2*k+1) * tcg * d6j                    
        end 
    end 
    Q5gh = sqrt(35.0/3.0) * (-1)^J * d3jinvnd * sum
    sum = 0.0
    for j = ell-1:2:ell+1
        if j <0;continue;end
        for k = j-1:2:j+1
            if k<0; continue;end
            tcg = (-1)^(ell+1-j) * cg1s[ell+1,j+1] * cg1s[k+1,j+1] * cg2s[k+1,ellp+1]
            d6j = d6_121[ell+1,j+1,k+1] * d6_222[k+1,ellp+1,ell+1]
            sum += y_fm2 * ndQW1s[ell+1][j+1][V_i,V_j] * (2*k+1) * tcg * d6j  
        end
    end 
    Q5gh += sqrt(35.0/3.0) * (-1)^J * d3jinvnd * sum
    ## term 5i & 5j
    sum = 0.0
    for jp = ellp-1:2:ellp+1
        if jp <0;continue;end
        for k=ell-1:2:ell+1
            if k <0;continue;end
            for j = k-1:2:k+1
                if j <0;continue;end
                tcg = (-1)^(jp+1-j) * cg1s[jp+1,j+1]*cg1s[k+1,j+1] 
                tcg *= cg1s[jp+1,ellp+1] * cg1s[k+1,ell+1]
                if tcg == 0.0;continue;end
                d69j = d6_121[jp+1,j+1,k+1] * d9_12[jp+1,ellp+1,k+1,ell+1]
                if d69j == 0.0;continue;end
                sum += XYT * ndQW1s[jp+1][j+1][V_j,V_i] * (2*jp+1) *(2*k+1) * tcg *d69j
            end 
        end 
    end 
    Q5ij = (-1)^J *sqrt(35.0/2.0) * d3jinvnd * sum
    sum = 0.0
    for jp = ell-1:2:ell+1
        if jp <0;continue;end
        for k=ellp-1:2:ellp+1
            if k <0;continue;end
            for j = k-1:2:k+1
                if j <0;continue;end
                tcg = (-1)^(jp+1-j) * cg1s[jp+1,j+1]*cg1s[k+1,j+1] 
                tcg *= cg1s[jp+1,ell+1] * cg1s[k+1,ellp+1]
                if tcg == 0.0;continue;end
                d69j = d6_121[jp+1,j+1,k+1] * d9_12[jp+1,ell+1,k+1,ellp+1]
                if d69j == 0.0;continue;end
                sum += XYT * ndQW1s[jp+1][j+1][V_i,V_j] * (2*jp+1) *(2*k+1) * tcg *d69j
            end 
        end 
    end 
    Q5ij += (-1)^J *sqrt(35.0/2.0) * d3jinvnd * sum
    Q5 = Q5a +Q5b + Q5c+Q5d + Q5ef + Q5gh + Q5ij 

    ## term6
    ## term 6a 
    Q6a = -1.0/3.0 * (x_fm2*QWL2_lp + y_fm2*QWL2_l -XYT * QWL2_J)
    ## term 6b
    tQ = ndQW2s[ellp+1][ell+1][V_j,V_i] *(x_fm2+y_fm2)
    tcg = cg2s[ellp+1,ell+1] *(-1)^ellp * sqrt(5.0) /hat(ell)
    Q6b = - 1.0/(3.0*sqrt(5.0)) * tQ * tcg * d3jinvnd 
    ## term 6c
    sum = 0.0
    for jp = ellp-1:2:ellp+1
        if jp <0;continue;end
        for j = ell-1:2:ell+1
            if j <0;continue;end
            tcg = (-1)^(jp+ell+ellp) * cg1s[ellp+1,jp+1] * cg1s[ell+1,j+1] * cg2s[jp+1,j+1]
            d6j = d6_21[jp+1,j+1,ell+1,ellp+1]
            sum += XYT * ndQW2s[jp+1][j+1][V_j,V_i] * tcg * d6j * hat(jp)
        end 
    end 
    Q6c = -1.0/(3.0*sqrt(15.0)) * d3jinvnd * sum 
    ## term 6d
    sum = 0.0
    for j = ellp-2:2:ellp+2
        if j < 0; continue;end
        tcg = (-1)^j *cg1s[1+1,2+1] * cg2s[j+1,ellp+1] *cg2s[j+1,ell+1] * d6_222[j+1,ellp+1,ell+1]
        sum += tcg *y_fm2 * (2*j+1) * ndQW2s[j+1][ell+1][V_j,V_i] *sqrt(5.0/(2*ell+1))
    end 
    Q6d = -sqrt(7.0/15.0) * d3jinvnd * sum

    sum = 0.0
    for j = ell-2:2:ell+2
        if j < 0; continue;end
        tcg = (-1)^j *cg1s[1+1,2+1] * cg2s[j+1,ellp+1] *cg2s[j+1,ell+1] * d6_222[j+1,ell+1,ellp+1]
        sum += tcg *x_fm2 * (2*j+1) * ndQW2s[ellp+1][j+1][V_j,V_i] *sqrt(5.0/(2*ellp+1))
    end 
    Q6e = -sqrt(7.0/15.0) * d3jinvnd * sum
    ## term 6f
    sum = 0.0
    for jp = ellp-1:2:ellp+1
        if jp <0;continue;end
        for j = ell-1:2:ell+1
            if j <0;continue;end
            tcg = (-1)^(j) * cg1s[jp+1,ellp+1] * cg1s[j+1,ell+1] * cg2s[jp+1,j+1]
            d9j = d9_12[jp+1,ellp+1,j+1,ell+1]
            sum += XYT * ndQW2s[jp+1][j+1][V_j,V_i] * tcg * d9j * (2*jp+1) * (2*j+1) /hat(j) *sqrt(5.0)
        end 
    end             
    Q6f = sqrt(7.0/15.0) *d3jinvnd * sum
    Q6 = Q6a + Q6b +Q6c +Q6d + Q6e + Q6f
    QT = (Q4TA + Q5 + Q6) * 2.0/3.0
    return QT
end

## V_c4: tensor diagonal part Eq.(B17+Erratum) 
function Vt_c4_d(V_i,V_j,ell,ellp,J,x_fm2,y_fm2,XYT,rho,mpi2,
                 F0X,F0Y,F1X,F1Y,F2X,F2Y,F3X,F3Y,
                 QL0,QL0p,QJ,QLm1,QLp1,
                 QWL0yx,QWL0_lp1,QWL0_lm1,QWL0_lp2yx,QWL0_lm2,
                 QWJ0p1,QWJ0m1,QWL2_l,QWL2_lp,QWL2_J,QW2_lm1,QW2_lp1,
                 QX1,QX1p,QX1J,QX1_m1,QX1_p1,
                 ndQW1s,ndQW2s,wsyms,d3jinv)
    cg1s = wsyms.cg1s;   cg2s = wsyms.cg2s
    d6_121 = wsyms.d6_121; d6_21  = wsyms.d6_21
    d6_222 = wsyms.d6_222; d9_12  = wsyms.d9_12
    Q4TA = (0.5 *rho -mpi2*(F0X+F0Y)) * (
        (x_fm2 +y_fm2)/XYT *QL0 
        -0.5 * (2*ell+3)/(2*ell+1) *QLm1 
        -0.5 * (2*ell-1)/(2*ell+1) *QLp1 )
    Q4TA += 0.5*XYT *(
        ((2*ell+1)^2 /((2*ell-1)*(2*ell+3)) -2.0) *QWL0yx 
        + (2*ell-1)*(ell+1)/((2*ell+1)*(2*ell+3)) *QWL0_lp2yx )
    if ell >= 2 
        Q4TA += 0.5*XYT *( (2*ell+3)*(ell-1)/((2*ell+1)*(2*ell-1)) *QWL0_lm2 )
    end 
    Q4TA += (F0Y-2.0*F1Y) *(
        +5.0/24.0 * XYT *delta(ell,1)
        +y_fm2/(2*XYT) *(-y_fm2+x_fm2+mpi2)*QL0
        -0.125/(2*ell+1) *(-y_fm2+x_fm2+mpi2)*((2*ell+3)*QLm1 +(2*ell-1)*QLp1) )
    Q4TA += (F0X-2.0*F1X) *(
        +5.0/24.0 * XYT *delta(ell,1)
        +x_fm2/(2*XYT) *(-x_fm2+y_fm2+mpi2)*QL0
        -0.125/(2*ell+1) *(-x_fm2+y_fm2+mpi2)*((2*ell+3)*QLm1 +(2*ell-1)*QLp1) )
    Q4TA += F3Y/3.0 *( 
        +5.0*XYT/8.0*delta(ell,1)
        +y_fm2/(2*XYT) * (-y_fm2+5.0*x_fm2+3.0*mpi2)*QL0 
        -0.125/(2*ell+1) *(y_fm2+3.0*x_fm2+3.0*mpi2)*((2*ell+3)*QLm1 +(2*ell-1)*QLp1))
    Q4TA += F3X/3.0 *( 
        +5.0*XYT/8.0*delta(ell,1)
        +x_fm2/(2*XYT) * (-x_fm2+5.0*y_fm2+3.0*mpi2)*QL0 
        -0.125/(2*ell+1) *(x_fm2+3.0*y_fm2+3.0*mpi2)*((2*ell+3)*QLm1 +(2*ell-1)*QLp1))    
    Q4TA += -(y_fm2*F2Y+x_fm2*F2X)/3.0 *( (x_fm2+y_fm2)/XYT *QL0 
        -((2*ell+3)*QLm1 +(2*ell-1)*QLp1)/(2*ell+1) )
    # term (5) in SY note
    Q5 = 0.0
    ## term 5a
    tQW1p = @views ndQW1s[ell+1][ell+2]
    Q5a = 1.0/3.0 *(y_fm2+x_fm2)*  ( 
        (2*ell-1)/(2*(2*ell+1)) * (tQW1p[V_i,V_j]+tQW1p[V_j,V_i])
    )
    if ell > 1 # not ell > 0!!
        tQW1n = @views ndQW1s[ell+1][ell]
        Q5a += 1.0/3.0 * (y_fm2+x_fm2) *(
            (2*ell+3)/(2*(2*ell+1)) * (tQW1n[V_i,V_j]+tQW1n[V_j,V_i]))
    end 
    Q5 += Q5a
    ## term 5b
    sum = 0.0
    for jp = ell-1:ell+1
        if jp <0;continue;end
        tQWs = @views ndQW1s[jp+1]
        jp1 = jp+1
        for j = ell-1:ell+1
            if j <0;continue;end
            j1 = j+1
            phase =  ifelse(j-1%2==0,1.0,-1.0)
            tQs = @views tQWs[j+1]
            for k=ell-1:ell+1
                if k<0 ;continue;end                                    
                tcg = phase * (cg1s[jp1,j1] * (cg1s[j+1,k+1] * (cg1s[ell+1,jp+1] * cg1s[ell+1,k+1])))
                if tcg == 0.0;continue;end
                twsym = d6_121[jp+1,j+1,k+1] * d6_21[jp+1,k+1,ell+1,ell+1]
                if twsym == 0.0;continue;end
                tQ = tQs[V_j,V_i] + tQs[V_i,V_j]
                fQ5b = twsym * tcg* hat(jp)*hat(j) 
                sum += fQ5b * 0.5*XYT * tQ 
            end 
        end 
    end                    
    Q5b = sqrt(2.0)/3.0  *d3jinv * sum 
    Q5 += Q5b
    ## term 5c
    Q5c = 1.0/3.0 * ( 
        (y_fm2+x_fm2) * QX1
        + 5.0 *XYT /12.0 * delta(ell,1) *(F0Y+F0X)
        - XYT/(2*(2*ell+1)) * ( (2*ell+3) *QX1_m1+ (2*ell-1) *QX1_p1 )  )
    Q5 += Q5c
    ## term 5d&5e
    sum = 0.0                       
    for j = ell-1:ell+1
        if j <0;continue;end
        for jp = j-1:j+1
            if jp <0 || jp > lmax+1;continue;end
            tcg = cg1s[jp+1,j+1] *cg1s[ell+1,j+1] * cg2s[jp+1,ell+1] 
            if tcg == 0.0;continue;end
            twsym = d6_121[jp+1,j+1,ell+1] * d6_222[jp+1,ell+1,ell+1]
            if twsym != 0.0
                tQWs = y_fm2 * ndQW1s[jp+1][j+1][V_j,V_i] + x_fm2 * ndQW1s[jp+1][j+1][V_i,V_j]
                sum += twsym * tcg* (2*jp+1) * tQWs * (-1.0)^(jp+1)
            end
            twsym = d6_121[ell+1,j+1,jp+1] * d6_222[jp+1,ell+1,ell+1]                                      
            if twsym != 0.0
                tQWs = x_fm2 * ndQW1s[ell+1][j+1][V_j,V_i] +y_fm2 * ndQW1s[ell+1][j+1][V_i,V_j]                                    
                sum += twsym * tcg * (2*jp+1) * tQWs * (-1.0)^(ell+1) 
            end
        end 
    end 
    Q5 += sqrt(35.0/3.0) *d3jinv * sum
    ## tern 5f
    sum = 0.0
    for jp = ell-1:ell+1
        if jp <0;continue;end
        for j = ell-1:ell+1
            if j <0;continue;end
            for k=ell-1:ell+1
                if k<0 ;continue;end                                    
                tcg1s = (-1.0)^(1-j) * cg1s[jp+1,j+1] *cg1s[j+1,k+1] * cg1s[jp+1,ell+1] * cg1s[k+1,ell+1]
                if tcg1s == 0.0;continue;end
                twsym = d6_121[jp+1,j+1,k+1] * d9_12[jp+1,ell+1,k+1,ell+1]
                if twsym == 0.0;continue;end
                tQWs = ndQW1s[jp+1][j+1][V_j,V_i] + ndQW1s[jp+1][j+1][V_i,V_j]
                fQ5f = twsym * tcg1s* (2*jp+1) * hat(k)*hat(j) 
                sum += fQ5f * XYT * tQWs 
            end 
        end 
    end                    
    Q5 += -sqrt(35.0/2.0) *d3jinv * sum

    ## term 6a        
    Q6 = 0.0
    Q6a = -1.0/3.0 *(2*(x_fm2+y_fm2) * QWL2_l -XYT/(2*(2*ell+1))*( (2*ell+3)*QW2_lm1 +(2*ell-1)*QW2_lp1))
    ## term 6b
    sum = 0.0
    for jp = ell-1:ell+1
        if jp <0 ;continue;end
        for j = ell-1:ell+1
            if j <0 ;continue;end                    
            tcg = cg1s[jp+1,ell+1] *cg1s[j+1,ell+1] * cg2s[jp+1,j+1] * sqrt(5.0)* (-1)^j /hat(j)
            if tcg == 0.0;continue;end
            twsym = d6_21[jp+1,j+1,ell+1,ell+1]
            if twsym == 0.0;continue;end
            prod = XYT * ndQW2s[jp+1][j+1][V_j,V_i]*twsym*tcg
            sum += (2*jp+1) * (2*j+1) * prod 
            ffac = -1.0/(15.0*sqrt(3.0)*(2*ell+1)*l2l[ell+1]) *(2*jp+1) * (2*j+1) 
        end 
    end 
    Q6b = -1.0/(15.0*sqrt(3.0)) / (2*ell+1) *d3jinv *sum  
    ## term 6c
    sum = 0.0                    
    for j = ell-2:2:ell+2
        if j < 0; continue;end
        tcg = cg2s[j+1,ell+1]^2 *sqrt(5.0)/hat(ell) * (-1)^ell
        d6j = d6_222[j+1,ell+1,ell+1]
        sum += (y_fm2*ndQW2s[j+1][ell+1][V_j,V_i]+
                x_fm2*ndQW2s[ell+1][j+1][V_j,V_i]) *tcg *d6j *(2*j+1)
    end
    Q6c = -sqrt(14.0/5.0)/3.0 *d3jinv * sum
    ## term 6d
    sum = 0.0                    
    for jp = ell-1:2:ell+1
        if jp < 0 ;continue;end
        for j = ell-1:ell+1
            if j < 0; continue;end
            tcg = cg1s[jp+1,ell+1] * cg1s[j+1,ell+1] *cg2s[jp+1,j+1] *sqrt(5.0)/hat(j) * (-1.0)^j
            d9j = d9_12[jp+1,ell+1,j+1,ell+1]
            sum += (XYT*ndQW2s[jp+1][j+1][V_j,V_i]) *tcg *d9j *(2*j+1) *(2*jp+1)
        end
    end
    Q6d = sqrt(7.0/15.0) *d3jinv * sum
    Q6 = Q6a + Q6b +Q6c + Q6d
    return  (Q4TA+Q5+Q6)*2.0/3.0 
end

## V_cD: central part  Eq.(B18)   
function Vc_cD(ell,XYT,rho,mpi2,ttit,ttis,F0X,F0Y,QL0)
    tmp = QL0*rho*mpi2/XYT 
    if ell == 0; tmp -= 2.0*mpi2*(F0X+F0Y); end
    tmp *= ttis * ttit /3.0
    if ell == 0; tmp += 3.0 * ( rho -2.0*mpi2 * (F0X+F0Y));end
    return tmp
end
## V_D: tensor term non-diagonal,  Eq.(B19)
## J=0 term is omitted, since it is not needed
function Vt_cD_nd(ell,ellp,x_fm2,y_fm2,XYT,rho,F0X,F0Y,F1X,F1Y,F3X,F3Y,
                  QL0,QL0p,QJ,tS12,ttit) 
    QDTA = 2.0/3.0 * (F0X-2.0*F1X+F3X +F0Y-2.0*F1Y+F3Y) *(
        x_fm2 *delta(ellp,0) + y_fm2 *delta(ell,0) ) 
    QDTB = -1.0/3.0 * ( (y_fm2 * QL0 + x_fm2 *QL0p)/XYT - QJ)
    return (QDTA + QDTB*rho) * tS12 * ttit
end
## V_D: tensor term diagonal,  Eq.(B20)
function Vt_cD_d(ell,x_fm2,y_fm2,XYT,rho,ttit,ttis,tS12,
                 F0X,F0Y,F1X,F1Y,F3X,F3Y,QL0,QL0p,QLm1,QLp1)
    QDTA = 2.0/3.0 * (F0X-2.0*F1X+F3X+F0Y-2.0*F1Y+F3Y) *(
        x_fm2 *delta(ell,0) + y_fm2 *delta(ell,0)
        -5.0/6.0 * XYT*delta(ell,1))
    QDTB = -1.0/3.0 * ( (y_fm2 * QL0 + x_fm2 *QL0p)/ XYT
                        -0.5 *( (2*ell+3)/(2*ell+1)* QLm1
                                + (2*ell-1)/(2*ell+1) *QLp1 ) )
    return (QDTA + QDTB*rho) * tS12 * ttit
end
    
function QWL(kF,k,kp,ell,ellp,p,ts,ws,mpi2,QLdict;is_xmul=true)
    dtdk3 = 0.5 *kF
    s = 0.0
    if ell < 0 || ellp < 0;return s;end
    deno = 2.0 * (2.0*pi)^2
    if p % 2 == 0 
        deno *= (k*kp)^div(p,2)
    elseif p==1 
        if is_xmul 
            deno *= kp  
        else    
            deno *= k
        end 
    else    
        println("QWL for p=$p is not supported now")
    end 
    for (i,t) in enumerate(ts)
        w = ws[i]
        k3 = 0.5 * (t+1.0)*kF 
        x = (k^2 + k3^2 + mpi2)/(2.0*k*k3)
        xp = (kp^2 + k3^2 + mpi2)/(2.0*kp*k3)
        tmp = dtdk3 * w * QL(x,ell,ts,ws,QLdict) * QL(xp,ellp,ts,ws,QLdict) 
        if p%2 == 0
            tmp *= k3^p
        elseif p==1
            tmp *= k3
            if is_xmul;tmp *= x;end
        end
        s += tmp
    end 
    return s /deno
end 

function QWls(kF,k,kp,ell,ellp,p,ts,ws,mpi2,QLdict)
    dtdk3 = 0.5 *kF
    s = 0.0
    deno = 2.0 * (2*pi)^2
    if p  == 0 
        deno *= k*kp
    else    
        println("QWls for p=$p is not supported now")
    end 
    if ell!=ellp; println("ell!=ellp in QWls!!");end
    for (i,t) in enumerate(ts)
        w = ws[i]
        k3 = 0.5 * (t+1.0)*kF 
        x = (k^2 + k3^2 + mpi2)/(2.0*k*k3)
        xp = (kp^2 + k3^2 + mpi2)/(2.0*kp*k3)
        s += dtdk3 * k3 * w * (k * QL(x,ell,ts,ws,QLdict) * (QL(x,ell-1,ts,ws,QLdict)-QL(x,ell+1,ts,ws,QLdict))
                               + kp* QL(xp,ell,ts,ws,QLdict) * (QL(xp,ell-1,ts,ws,QLdict)-QL(xp,ell+1,ts,ws,QLdict)))
    end 
    return s/deno
end 

function QX(l,ndQW1s,i,j)
    s = 0.0
    if l !=0
        s  = l * ndQW1s[l+1][l][j,i] 
    end 
    if l+2 <= lmax + 1
        s += (l+1) * ndQW1s[l+1][l+2][j,i]
    end 
    return s / (2.0*l+1.0)
end 
