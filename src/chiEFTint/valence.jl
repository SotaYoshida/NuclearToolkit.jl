# pfunc(p,p',P)
# 1: P^2 2: P^2 3: (ddp +ddp') 4: dd p-dd p' 5: P^2
function calc_pPfac(xr,wr,xrP,wrP,n,ell,np,ellp,N,L,Np,Lp,Rnl,RNL,pnrank,dwn,pfunc,to)
    Rnlx = @views Rnl[n+1,ell+1,:]
    Rnly = @views Rnl[np+1,ellp+1,:]
    RNL1 = @views RNL[N+1,L+1,:]
    RNL2 = @views RNL[Np+1,Lp+1,:]    
    r = 0.0
    @inbounds for (i,x) in enumerate(xr)
        xdwn2= (x*dwn)^2
        ex = sqrt(1.0+xdwn2)
        xRw = Rnlx[i]* wr[i] * xdwn2
        @inbounds for (j,y) in enumerate(xr)
            ydwn2= (y*dwn)^2
            ey = sqrt(1.0+ydwn2)
            ree = 1.0/sqrt(ex*ey)
            yRw = Rnly[j]* wr[j] *ydwn2
            @inbounds for (k,P) in enumerate(xrP)
                tRNL = RNL1[k] * RNL2[k] * (P *dwn)^2 * wrP[k] 
                PRw = pfunc(x,ell,y,ellp,P) * tRNL                 
                r += xRw * yRw* ree * PRw
            end
        end
    end
    return r 
end

function NLOvs(chiEFTobj,dLECs,vs,xr,wr,xrP,wrP,Rnl,RNL,
               cg1s,sp_P5_9j,nljsnt,pnrank,ip,X9,U6,t5v,
               Jtot,iza,ia,izb,ib,izc,ic,izd,id,to)
    emax = chiEFTobj.emax
    U6_j = U6[Jtot+1]
    na,la,jda = nljsnt[ia]; nb,lb,jdb = nljsnt[ib]
    nc,lc,jdc = nljsnt[ic]; nd,ld,jdd = nljsnt[id]
    mab = 2*na+la+2*nb+lb; mcd = 2*nc+lc+2*nd+ld
    Nab = Nfac_jj(na,la,jda,nb,lb,jdb)
    Ncd = Nfac_jj(nc,lc,jdc,nd,ld,jdd)
    lrmax = jmax + 1
    elab=2*na+la+2*nb+lb; elabp=2*nc+lc+2*nd+ld
    Nrmax = 2*emax
    MN = Ms[pnrank];dwn = 1.0/MN
    mfac = 1.0    
    rc_vs_1 = dLECs["c_vs_1"] * mfac
    rc_vs_2 = dLECs["c_vs_2"] * mfac
    rc_vs_3 = dLECs["c_vs_3"] * mfac
    rc_vs_4 = dLECs["c_vs_4"] * mfac
    rc_vs_5 = dLECs["c_vs_5"] * mfac
    for S=0:1
        tX9 = X9[S+1]
        U6_S = U6_j[S+1]
        for Sp=0:1
            U6_Sp = U6_j[Sp+1]
            tX9p = X9[Sp+1]
            if S ==Sp # term 1,2,5
                ell = ellp = 0
                pfunc = fp_P2
                minLam = max(abs(Jtot-S),abs(la-lb))
                maxLam = min(Jtot+S,la+lb)
                if minLam > maxLam;continue;end
                for Lam = minLam:maxLam
                    L = Lam                        
                    minLamp = max(abs(Jtot-Sp),abs(lc-ld))
                    maxLamp = min(Jtot+Sp,lc+ld)
                    if minLamp > maxLamp;continue;end
                    for Lamp = minLamp:maxLamp
                        Lp = Lamp
                        nja = jda-2*la; njb = jdb-2*lb
                        t5v[1]=la;t5v[2]=nja;t5v[3]=lb
                        t5v[4]=njb;t5v[5]=Lam
                        x1 = get(tX9[Jtot+1],t5v,0.0)
                        
                        njc = jdc-2*lc; njd = jdd-2*ld
                        t5v[1]=lc;t5v[2]=njc;t5v[3]=ld
                        t5v[4]=njd;t5v[5]=Lamp
                        x2 = get(tX9p[Jtot+1],t5v,0.0)
                        if x1==0.0 || x2 == 0.0;continue;end
                        t9j = x1*x2* ifelse((Lam+Lamp)%2==0,1.0,-1.0)
                        t6j = 1.0 # due to ell=ellp=0
                        for N=0:div(elab,2)
                            if (elab-L) % 2 !=0;continue;end
                            n = div(elab-(2*N+L),2);if n<0;continue;end
                            for Np=0:div(elabp,2)
                                if (elabp-Lp) % 2 !=0;continue;end
                                np = div(elabp-(2*Np+Lp),2); if np<0;continue;end
                                HOB1=HObracket(na,la,nb,lb,n,ell,N,L,Lam)
                                HOB2=HObracket(nc,lc,nd,ld,np,ellp,Np,Lp,Lamp)
                                wsyms = t9j * t6j
                                coeff = HOB1*HOB2 * wsyms
                                if coeff == 0.0;continue;end
                                pPfac = calc_pPfac(xr,wr,xrP,wrP,n,ell,np,ellp,N,L,Np,Lp,Rnl,RNL,pnrank,dwn,pfunc,to)
                                for Jrel =0:jmax
                                    Jrelp = Jrel
                                    if Jrel == S
                                        rv1 = rc_vs_1 * coeff * 4*pi * pPfac
                                        rv2 = (2*S*(S+1)-3) * rc_vs_2 * coeff * 4*pi * pPfac 
                                        vs[1] += rv1
                                        vs[2] += rv2
                                    end                                    
                                    sum5 = 0.0
                                    for a=0:2:2
                                        t = hat(a)
                                        t *= cg1s[div(a,2)+1]
                                        t *= clebschgordan(Float64,Lp,0,a,0,L,0)
                                        t *= wigner6j(Float64,Jrel,L,Jtot,Lp,Jrelp,a)
                                        t *= sp_P5_9j[div(a,2)+1,S+1]
                                        t *= (-1.0)^(L+Jtot+Jrel)
                                        sum5 += t
                                    end
                                    fac5 = 24*pi *hat(Lp) * hat(S)^2 * sum5
                                    vs[5] += rc_vs_5 * fac5 * pPfac                                    
                                end
                            end
                        end
                    end
                end
            else
                ## term 3,4 S+Sp=1
                pfunc = fp_P2
                for ellidx=1:2
                    ell = ellp = 0
                    if ellidx == 1
                        ell = 1 ; ellp = 0 # p' ~0
                    else
                        ell = 0 ; ellp = 1 # p ~ 0
                    end
                    minLam = max(abs(Jtot-S),abs(la-lb))
                    maxLam = min(Jtot+S,la+lb)
                    if minLam > maxLam;continue;end
                    for Lam = minLam:maxLam
                        minLamp = max(abs(Jtot-Sp),abs(lc-ld))
                        maxLamp = min(Jtot+Sp,lc+ld)
                        if minLamp > maxLamp;continue;end
                        for L = abs(Lam-ell):Lam+ell
                            for Lamp = minLamp:maxLamp
                                for Lp = abs(Lamp-ellp):Lamp+ellp
                                    nja = jda-2*la; njb = jdb-2*lb
                                    njc = jdc-2*lc; njd = jdd-2*ld
                                    t5v[1]=la;t5v[2]=nja;t5v[3]=lb;
                                    t5v[4]=njb;t5v[5]=Lam
                                    x1 = get(tX9[Jtot+1],t5v,0.0)
                                    t5v[1]=lc;t5v[2]=njc;t5v[3]=ld;
                                    t5v[4]=njd;t5v[5]=Lamp
                                    x2 = get(tX9p[Jtot+1],t5v,0.0)
                                    if x1*x2 == 0.0;continue;end
                                    t9j= x1*x2*ifelse((Lam+Lamp)%2==0,1.0,-1.0)
                                    for N=0:div(elab,2)
                                        if (elab-L) % 2 !=0;continue;end
                                        n = div(elab-(2*N+L),2);if n<0;continue;end
                                        for Np=0:div(elabp,2)
                                            if (elabp-Lp) % 2 !=0;continue;end
                                            np = div(elabp-(2*Np+Lp),2); if np<0;continue;end
                                            HOB1=HObracket(na,la,nb,lb,n,ell,N,L,Lam)
                                            HOB2=HObracket(nc,lc,nd,ld,np,ellp,Np,Lp,Lamp)
                                            if HOB1 * HOB2 == 0.0; continue;end
                                            #@timeit to "pP"
                                            pPfac = calc_pPfac(xr,wr,xrP,wrP,n,ell,np,ellp,N,L,Np,Lp,Rnl,RNL,pnrank,dwn,pfunc,to)
                                            for Jrel = abs(ell-S):ell+S
                                                for Jrelp = abs(ellp-Sp):ellp+Sp
                                                    u6j =U6_S[Lam+1][L+1][ell+1][Jrel-abs(ell-S)+1]
                                                    u6jp=U6_Sp[Lamp+1][Lp+1][ellp+1][Jrelp-abs(ellp-Sp)+1]
                                                    t6j = u6j * u6jp
                                                    coeff = HOB1*HOB2 * t6j * t9j
                                                    if coeff == 0.0;continue;end                                                    
                                                    fac34 = 12.0*pi *hat(Jrel) * hat(Jrelp) * hat(Lp)
                                                    fac34 *= clebschgordan(Float64,Lp,0,1,0,L,0)
                                                    fac34 *= wigner6j(Float64,Jrel,L,Jtot,Lp,Jrelp,1)
                                                    fac34 *= wigner9j(Float64,1,1,1,ell,S,Jrel,ellp,Sp,Jrelp)
                                                    fac4 = fac34 
                                                    fac3 = 2.0*sqrt(2.0) * fac34
                                                    fac3 *= (-1)^(L+Jtot+Jrelp+Sp+1)
                                                    fac4 *= (-1)^(L+Jtot+Jrelp)
                                                    vs[3] += rc_vs_3 * fac3 * pPfac                                    
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
    if pnrank != 2
        vs .*= Nab * Ncd
    end
    return nothing
end

