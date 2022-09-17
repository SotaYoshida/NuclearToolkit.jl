"""
    OPEP(chiEFTobj,to;pigamma=true,debugmode=false)

calc. One-pion exchange potential in the momentum-space

Reference: R. Machleidt, Phys. Rev. C 63 024001 (2001).
"""
function OPEP(chiEFTobj,to;pigamma=true,debugmode=false)
    ts = chiEFTobj.ts; ws = chiEFTobj.ws; xr = chiEFTobj.xr; V12mom = chiEFTobj.V12mom
    dict_pwch = chiEFTobj.dict_pwch; 
    lsjs = chiEFTobj.lsjs; tllsj = chiEFTobj.tllsj
    opfs = chiEFTobj.opfs; QLdict = chiEFTobj.util_2n3n.QWs.QLdict
    n_mesh = chiEFTobj.params.n_mesh
    hc3 = hc^3
    tVs = zeros(Float64,6)
    opfs = zeros(Float64,8)
    mpi0 = mpis[2]; mpi02 = mpi0^2
    mpipm = mpis[1]; mpipm2 = mpipm^2
    for pnrank =1:3
        tdict = dict_pwch[pnrank]
        MN = Ms[pnrank];dwn = 1.0/MN;sq_dwn=dwn^2
        fff = pi / ((2*pi)^3 * MN^2)
        coeff = -(MN*gA/(2*Fpi))^2
        itt = itts[pnrank]
        tllsj[1] = itt
        @inbounds for J=0:jmax
            lsj = lsjs[J+1]
            f_idx = 6
            if J==0; f_idx = 3;end                
            @inbounds for i= 1:n_mesh
                x = xr[i];xdwn = x * dwn;sq_xdwn= xdwn^2
                ex = sqrt(1.0+sq_xdwn)
                @inbounds for j = 1:n_mesh
                    y = xr[j]; ydwn = y*dwn;sq_ydwn= ydwn^2
                    ey = sqrt(1.0+sq_ydwn)
                    nfac = 1.0/(x* y* sq_dwn)
                    ree = 1.0/sqrt(ex*ey) * freg(x,y,4)
                    f_sq!(opfs,xdwn,ydwn)
                    if pnrank != 2 ## pp/nn
                        cib_lsj_opep(opfs,x,y,mpi02,1,J,pnrank,nfac,ts,ws,tVs,QLdict)
                    else ## pn
                        cib_lsj_opep(opfs,x,y,mpi02,1,J,pnrank,nfac,ts,ws,tVs,QLdict)
                        cib_lsj_opep(opfs,x,y,mpipm2,2,J,pnrank,nfac,ts,ws,tVs,QLdict;additive=true)
                        #if pigamma
                        #    cib_lsj_opep(opfs,x,y,mpipm2,2,J,pnrank,nfac,ts,ws,tVs,pigamma;additive=true)
                        #end
                    end
                    t_fc = hc3 * fff * coeff * ree
                    @inbounds for idx = 1:f_idx
                        @views tllsj[2:5] .= lsj[idx]
                        tl,tlp,tS,tJ = lsj[idx] 
                        if pnrank%2 == 1 && (tl+tS+1)%2 != 1;continue;end
                        V12idx = get(tdict,tllsj,-1)
                        if V12idx == -1;continue;end                        
                        tfac = tVs[idx] * t_fc
                        V12mom[V12idx][i,j] += tfac
                    end
                end
            end
        end
    end
    return nothing
end

function fac_pig(beta,c5=0.0)
    return - (1.0-beta)^2 / (2*beta^2) * log(1+beta) +(1.0+beta)/(2*beta) -2.0*c5
end

function cib_lsj_opep(opfs,x,y,mpi2,nterm,J,pnrank,facin,ts,ws,tVs,QLdict,pigamma=false;additive=false)
    x2 = x^2; y2 = y^2
    z = (mpi2+x2+y2) / (2.0*x*y)
    QJ = QJm1 = 0.0
    nfac = facin
    if pigamma
        nfac = facin * fsalpha/pi
        q2s = zeros(Float64,length(ts))
        for (i,t) in enumerate(ts)
            q2 = x2 + y2 -2.0*x*y*t
            beta = q2/mpi2
            q2s[i] = fac_pig(beta) * t
        end
        QJ = QL(z,J,q2s,ws,QLdict) 
        if J>0;QJm1=QL(z,J-1,q2s,ws,QLdict);end
    else
        QJ = QL(z,J,ts,ws,QLdict)
        if J>0;QJm1=QL(z,J-1,ts,ws,QLdict);end
    end
    IJ0 = nfac * QJ
    IJ1 = nfac * (z * QJ -delta(J,0)) #Eq. (B19)
    IJ2 = nfac * (J*z* QJ + QJm1) /(J+1) #Eq. (B20)
    IJ3 = nfac * sqrt(J/(J+1)) * (z* QJ - QJm1) #Eq. (B21)
    #Eq. (B28)
    v1 = opfs[1] * IJ0 + opfs[2] *IJ1 
    v2 = opfs[3] * IJ0 + opfs[4] *IJ2
    v3 = opfs[5] * IJ0 + opfs[6] *IJ1
    v4 = opfs[4] * IJ0 + opfs[3] *IJ2
    v5 = opfs[7] * IJ3
    v6 = -v5                        
    if J==0; v2=v4=v5=v6=0.0;end
    if J%2==1 && pnrank!=2; v1 =0.0;end
    if J%2==0 && pnrank!=2; v2 =0.0;end
    if J%2!=0 && pnrank!=2; v3=v4=v5=v6=0.0;end
    v34 = -sqrt(J*(J+1)) *(v3-v4)
    v56 = sqrt(J*(J+1)) * (v5+v6)
    d2j1 = 1.0/(2*J+1)
    if nterm == 1
        phase = ifelse(pnrank==2,-1.0,1.0)         
        tVs[1] = additive_sum(additive,tVs[1],v1 *phase)
        tVs[2] = additive_sum(additive,tVs[2],v2 *phase)
        tVs[3] = additive_sum(additive,tVs[3],d2j1 * ((J+1)* v3 + J*v4-v56)*phase)
        tVs[4] = additive_sum(additive,tVs[4],d2j1 * ( J*v3 + (J+1)*v4 +v56) *phase)
        tVs[5] = additive_sum(additive,tVs[5],-d2j1 * (v34-(J+1)*v5+J*v6)*phase)
        tVs[6] = additive_sum(additive,tVs[6],-d2j1 * (v34+J*v5-(J+1)*v6)*phase)
    else
        is = J%2 + 1
        it = is%2 +1
        ttis = ifelse(is==2,-2.0,2.0)
        ttit = ifelse(it==2,-2.0,2.0)
        tVs[1] = additive_sum(additive,tVs[1],ttis * v1)
        tVs[2] = additive_sum(additive,tVs[2],ttit * v2)
        tVs[3] = additive_sum(additive,tVs[3],d2j1 * ((J+1)* (ttis*v3) + J*(ttis*v4)-(ttis*v56)))
        tVs[4] = additive_sum(additive,tVs[4],d2j1 * ( J*(v3*ttis) + (J+1)*(ttis*v4) +(ttis*v56)))
        tVs[5] = additive_sum(additive,tVs[5],-d2j1 * ((ttis*v34)-(J+1)*(ttis*v5)+J*(ttis*v6)))
        tVs[6] = additive_sum(additive,tVs[6],-d2j1 * ((ttis*v34)+J*(ttis*v5)-(J+1)*(ttis*v6)))
    end
    return nothing 
end

function additive_sum(TF::Bool,retv,inv)
    if TF
         return retv + inv 
    else
        return inv
    end
end

### function-forms for partial waves
f1(x,y,c,cdum) = c
fx(x,y,c,cdum) = c * x
fy(x,y,c,cdum) = c * y
fxxpyy(x,y,c,cdum) = c * (x^2 + y^2)
fxy(x,y,c,cdum) = c *x*y
fx2(x,y,c,cdum) = c * x^2
fy2(x,y,c,cdum) = c * y^2
f_442(x,y,c1,c2) = c1 * (x^4 + y^4) + c2 *(x^2 * y^2)
f_x42(x,y,c1,c2) = c1 * x^4  + c2 *(x^2 * y^2)
f_y42(x,y,c1,c2) = c1 * y^4  + c2 *(x^2 * y^2)
f_x2y2(x,y,c,cdum) = c * x^2 * y^2
f_x3y(x,y,c,cdum) = c * x^3 * y
f_x3y3(x,y,c,cdum) = c * x^3 * y^3
f_xy3(x,y,c,cdum) = c * x * y^3
f_31(x,y,c,cdum) = c * (x^3 * y + x * y^3)
fp_P2(p,ell,pp,ellp,P) = P^2
fp_ddP(p,ell,pp,ellp,P) = P * (delta(ell,1)*delta(ellp,0)*p + delta(ell,0)*delta(ellp,1)*pp)

function set_pjs!(J,pjs,ts)
    pjs .= 0.0
    if J ==0
        pjs[:,1] .= 1.0; pjs[:,3] .= 0.0
    else
        pjs[:,3] .= 1.0
        for (i,t) in enumerate(ts)        
            pjs[i,1] = t
        end
        if J>1
            for (i,t) in enumerate(ts)
                pj = pjs[i,1]
                pjm1 = pjs[i,3]            
                for tJ = 2:J
                    a = t * pj
                    b = a-pjm1
                    pjm1 = pj
                    pj = -b/tJ + b+a
                end
                pjs[i,1] = pj
                pjs[i,3] = pjm1
            end
        end
    end
    for (i,t) in enumerate(ts)
        pjs[i,2] = pjs[i,1] * t
        pjs[i,4] = pjs[i,2] * t
        pjs[i,6] = pjs[i,4] * t
        pjs[i,5] = pjs[i,3] * t
        pjs[i,7] = pjs[i,5] * t
    end

    return nothing
end

"""

`9`x`nthreads` matrix to store sum of each tpe channal, C/T/S/LS/SigmaL
"""
mutable struct tpe_ch
    Vc::Vector{Float64}
    Wc::Vector{Float64}
    Vt::Vector{Float64}
    Wt::Vector{Float64}
    Vs::Vector{Float64}
    Ws::Vector{Float64}
    Vls::Vector{Float64}
    Wls::Vector{Float64}
    Vsl::Vector{Float64}
end
"""
    tpe(chiEFTobj,to::TimerOutput)

calc. two-pion exchange terms up to N3LO(EM) or N4LO(EMN)

The power conting schemes for EM/EMN are different;
The ``1/M_N`` correction terms appear at NNLO in EM and at N4LO in EMN.

# References
- EM: R. Machleidt and D.R. Entem [Physics Reports 503 (2011) 1–7](https://doi.org/10.1016/j.physrep.2011.02.001)
- EMKN: D. R. Entem, N. Kaiser, R. Machleidt, and Y. Nosyk, [Phys. Rev. C 91, 014002 (2015)](https://doi.org/10.1103/PhysRevC.91.014002).
"""
function tpe(chiEFTobj,to)
    LECs = chiEFTobj.LECs.dLECs
    ts = chiEFTobj.ts
    dict_pwch = chiEFTobj.dict_pwch
    lsjs = chiEFTobj.lsjs; tllsj = chiEFTobj.tllsj; opfs = chiEFTobj.opfs
    nthre = nthreads()
    c1_NNLO = LECs["c1_NNLO"];c2_NNLO = LECs["c2_NNLO"];c3_NNLO = LECs["c3_NNLO"];c4_NNLO = LECs["c4_NNLO"]
    d12 = LECs["d12"];d3 = LECs["d3"]; d5 = LECs["d5"]; d145 = LECs["d145"];e14 = LECs["e14"];e17 = LECs["e17"]
    mmpi = sum(mpis)/3.0
    pjs = zeros(Float64,length(ts),7)
    pjs_para = [ deepcopy(pjs) for i=1:nthre]
    
    tmpsum = [ [zeros(Float64,7) for j=1:9] for i=1:nthreads()]
    tmpLECs = Dict{String,Float64}()
    tmpLECs["c1"] = tmpLECs["c2"] = tmpLECs["c3"] = tmpLECs["c4"] = 0.0 
    tmpLECs["r_d12"] = tmpLECs["r_d3"] = tmpLECs["r_d5"] = tmpLECs["r_d145"] = 0.0
    opfs_para = [ deepcopy(opfs) for i=1:nthre] 
    gis_para = [ [zeros(Float64,7) for i=1:9] for j=1:nthre]
    tVs_para = [ zeros(Float64,6) for i=1:nthre]   

    for pnrank =1:3
        tdict = dict_pwch[pnrank]
        MN = Ms[pnrank];dwn = 1.0/MN
        fff = pi / ((2*pi)^3 * MN^2)
        nd_mpi = mmpi/MN        
        Fpi2 = (Fpi/MN)^2
        c1 = c1_NNLO * MN * 1.e-3
        c2 = c2_NNLO * MN * 1.e-3
        c3 = c3_NNLO * MN * 1.e-3
        c4 = c4_NNLO * MN * 1.e-3
        r_d12 = d12 * MN^2 * 1.e-6
        r_d3 = d3 * MN^2 * 1.e-6
        r_d5 = d5 * MN^2 * 1.e-6
        r_d145 = d145 * MN^2 * 1.e-6
        r_e14 = e14 * MN^3 * 1.e-9 
        r_e17 = e17 * MN^3 * 1.e-9
        tmpLECs["c1"] = c1
        tmpLECs["c2"] = c2
        tmpLECs["c3"] = c3
        tmpLECs["c4"] = c4 
        tmpLECs["r_d12"] = r_d12
        tmpLECs["r_d3"] = r_d3
        tmpLECs["r_d5"] = r_d5
        tmpLECs["r_d145"] = r_d145
        tmpLECs["r_e14"] = r_e14
        tmpLECs["r_e17"] = r_e17
        itt = itts[pnrank]
        tllsj[1] = itt
        tllsj_para = [ deepcopy(tllsj) for i=1:nthre]
        LamSFR_nd = chiEFTobj.params.LambdaSFR * dwn
        LoopObjects = precalc_2loop_integrals(chiEFTobj,LamSFR_nd,nd_mpi,Fpi2,c1,c2,c3,c4,r_d12,r_d3,r_d5,r_d145,r_e14,r_e17)
        @inbounds for J=0:jmax
            lsj = lsjs[J+1]
            f_idx = 6
            if J==0; f_idx = 3;end
            set_pjs!(J,pjs,ts)
            for i=1:nthre; pjs_para[i] .= pjs;end
            tpe_for_givenJT(chiEFTobj,LoopObjects,Fpi2,tmpLECs,
                            J,pnrank,fff,dwn,nd_mpi,pjs_para,gis_para,opfs_para,
                            f_idx,tVs_para,lsj,tllsj_para,tdict,tmpsum,to)
        end
    end
    return nothing
end

function f_sq!(opf,xdwn,ydwn)
    opf .= 0.0
    opf[1] = opf[6] = -2.0 * (xdwn^2 + ydwn^2)
    opf[2] = opf[5] = 4.0 * xdwn * ydwn
    opf[3] = -opf[1]; opf[4] = -opf[2]
    opf[7] = 2.0*( xdwn^2 - ydwn^2)
    opf[8] = -opf[7]
    return nothing
end

function f_ss!(opf)
    opf .= 0.0
    opf[1] = -6.0
    opf[2] = opf[4] = opf[5] = 0.0
    opf[3] = opf[6] = 2.0
    opf[7] = opf[8] = -2.0
    return nothing
end

function f_c!(opf)
    opf .= 0.0
    opf[1] = opf[3] = opf[6] = 2.0
    opf[2] = opf[4] = opf[5] = 0.0
    opf[7] = opf[8] = -2.0
    return nothing
end

function f_ls!(opf,xdwn,ydwn)
    xy2 = 2.0 * xdwn * ydwn
    opf .= 0.0
    opf[4] = opf[5] = opf[11] = -xy2
    opf[10] = xy2
    return nothing
end

function f_sl!(opf,xdwn,ydwn)
    xxyy2 = 2.0 * xdwn^2 * ydwn^2
    opf .= 0.0
    opf[1] = opf[3] = opf[7] = opf[8] = -xxyy2
    opf[6] = xxyy2
    opf[9] = 2.0 * xxyy2
    return nothing
end

"""
    calc_LqAq(w,q,nd_mpi,usingSFR,LamSFR)
To calculate L(q) and A(q) appear in TPE contributions; L(q)=EM Eq.(4.11), A(q)=EM Eq.(D.28).
If usingSFR, Lq&Aq are replaced by Eqs.(2.24)&(C3) in EKMN paper.
"""
function calc_LqAq(w,q,nd_mpi,usingSFR,LamSFR)
    if usingSFR
        s = sqrt(LamSFR^2-4*nd_mpi^2)
        Lq = w/(2*q) * log((LamSFR^2 * w^2 + q^2 * s^2 + 2*LamSFR*q*w*s)/(4.0*nd_mpi^2 * (LamSFR^2 + q^2) ) )
        Aq = atan(q*(LamSFR-2*nd_mpi)/(q^2 + 2*LamSFR*nd_mpi))  / (2.0*q)  
        return Lq,Aq
    else
        Lq = w/q * log((w+q)/(2.0*nd_mpi) )
        Aq = atan(q/(2.0*nd_mpi))  / (2.0*q)
        return Lq,Aq
    end
end

"""
    Vt_term(chi_order,w,q2,k2,Lq,Aq,nd_mpi,r_d145,Fpi4)

- NLO: EM Eq.(4.10)
- NNLO: EM Eq.(4.15) & Eq.(4.22) => `it_pi` term
- N3LO: EM Eq.(D.11) (`1/M^2_N`` term), Eq.(D.23) (2-loop term)
"""
function Vt_term(chi_order,LoopObjects,w,tw2,q2,k2,Lq,Aq,nd_mpi,r_d145,Fpi2;EMN=false,calc_TPE_separately=false)
    Fpi4 = Fpi2^2; w2 = w^2; tw2Aq = tw2 * Aq
    tmp_s = 0.0
    # NLO
    if chi_order >= 1
        f_NLO_Vt = -3.0* gA4 / (64.0 * pi^2 * Fpi4)       
        tmp_s += Lq * f_NLO_Vt
    end
    # NNLO
    if chi_order >= 2
        f_NNLO_Vt = 3.0 * gA4 / (256.0 * pi * Fpi4)
        if EMN
            if chi_order >= 3
                tmp_s += (5*nd_mpi^2 + 2* q2) * Aq * f_NNLO_Vt
            end
        else
            it_pi = (nd_mpi + w2 * Aq) ##nd_mpi term is different from EMN
            tmp_s += (3*tw2Aq +it_pi) * f_NNLO_Vt/2 
        end
    end
    # N3LO
    if chi_order >= 3 
        if !EMN
            ## EM eq.(D.11)
            f_N3LO_Vt = gA4 / (32.0 * pi^2 * Fpi4) 
            tmp_s += Lq * (k2 +3.0/8.0 *q2 + nd_mpi^4 /w2) *f_N3LO_Vt
            ## EM eq.(D.23)
            f_N3LO_2l_Vt = -gA2 * r_d145 /(32.0*pi^2 *Fpi4)
            tmp_s += w2 * Lq * f_N3LO_2l_Vt
        else           
            obj = LoopObjects.n3lo; tmp_s += Integral_ImV_T(obj.mudomain,q2,obj,obj.ImVt)
        end
    end
    if chi_order >= 4 && EMN && calc_TPE_separately
        obj = LoopObjects.n4lo; tmp_s += Integral_ImV_T(obj.mudomain,q2,obj,obj.ImVt)  # (2.15) VtB
        obj = LoopObjects.n4lo; tmp_s += Integral_ImV_T(obj.mudomain3,q2,obj,obj.ImVt3) #  3loop Vt12(2.27)&Vt13(2.32) 
    end
    return tmp_s
end

function Wt_term(chi_order,LoopObjects,w,q2,Lq,Aq,nd_mpi,c4,Fpi2;EMN=false,calc_TPE_separately=false)
    Fpi4 = Fpi2^2; Fpi6 = Fpi2^3; w2 = w^2; w2Aq = w2 * Aq
    nd_mpi2 = nd_mpi^2; nd_mpi4 = nd_mpi2^2
    tmp_s = 0.0
    if chi_order >=2
        #NNLO
        fac_NNLO = - gA2 / (32.0 * pi * Fpi4)
        if EMN
            ### EMN version non 1/M
            tmp_s += (c4 * w2 *Aq ) * fac_NNLO 
            if chi_order >=3
                tmp_s += -1/4 * (gA2*(3*nd_mpi2 + q2) -w2) *Aq *fac_NNLO
            end
        else
            ## EM eq.(4.15), nd_mpi in it-2PE term is different from EKMN
            term1 = (c4+0.25) * w2
            term2 = -0.125*gA2 * (10.0 *nd_mpi2 + 3.0 *q2)
            it_pi = gA2 /8.0 * (nd_mpi + w2Aq)
            tmp_s += ((term1+term2)*Aq +it_pi) * fac_NNLO            
        end
    end
    #N3LO        
    if chi_order >= 3        
        f_N3LO_a_Wt = c4^2 / (96.0*pi^2 *Fpi4)
        tmp_s += Lq * w2 * f_N3LO_a_Wt       
        f_N3LO_b_Wt = - c4 / (192.0*pi^2 *Fpi4)
        tmp_s += Lq * ( gA2 *( 16.0 * nd_mpi2 + 7.0 * q2) -w2) * f_N3LO_b_Wt # ci/MN
        f_N3LO_c_Wt = 1.0/ (1536.0* pi^2 * Fpi4) 
        f_N3LO_2l_Wt= gA4 / (2048.0 * pi^2 * Fpi6)    
        if EMN
            obj = LoopObjects.n3lo; tmp_s += Integral_ImV_T(obj.mudomain,q2,obj,obj.ImWt) #2loop 
        else
            tmp_s += Lq * ( 4.0*gA4 *(7.0*nd_mpi2 + 17.0/4.0 * q2 + 4.0 * nd_mpi4 /w2)
                           -32.0*gA2 *(nd_mpi2 + 7.0/16.0 * q2) + w2) * f_N3LO_c_Wt
            tmp_s += f_N3LO_2l_Wt * w2Aq *( w2Aq + 2.0*nd_mpi * (1.0+2.0*gA2))
        end
    end
    if chi_order >= 4 && EMN && calc_TPE_separately
        obj = LoopObjects.n4lo; tmp_s += Integral_ImV_T(obj.mudomain,q2,obj,obj.ImWt)
        obj = LoopObjects.n4lo; tmp_s += Integral_ImV_T(obj.mudomain3,q2,obj,obj.ImWt3)# 3loop Wt13 EKMN (2.29) Wt13 EKMN (2.34)
    end
    return tmp_s
end

function Vs_term(chi_order,LoopObjects,w,tw2,q2,k2,Lq,Aq,nd_mpi,r_d145,Fpi2;EMN=false,calc_TPE_separately=false)
    Fpi4 = Fpi2^2; nd_mpi4 = nd_mpi^4; w2 = w^2; tw2Aq = tw2 * Aq 
    tmp_s = 0.0
    ## NLO
    if chi_order >=1
        f_NLO_Vs = -3.0 * gA4  / (64.0 * pi^2 * Fpi4)
        tmp_s += Lq * f_NLO_Vs 
    end
    ## NNLO        
    if chi_order >= 2
        f_NNLO_Vs = 3.0 * gA4 / (256.0 * pi * Fpi4)
        if EMN 
            if chi_order >= 3
                tmp_s += (5*nd_mpi^2 + 2* q2) * Aq * f_NNLO_Vs
            end
        else #EM
            it_pi = (nd_mpi + w2 * Aq)  ##nd_mpi term is different from EKMN
            tmp_s += ( 3 * tw2Aq +it_pi) * f_NNLO_Vs/2 
        end
    end
    ## N3LO 
    if chi_order >= 3 
        if EMN
            obj = LoopObjects.n3lo; tmp_s += Integral_ImV_T(obj.mudomain,q2,obj,obj.ImVt) 
        else
            f_N3LO_Vs = gA4 / (32.0 * pi^2 * Fpi4)
            f_N3LO_2l_Vs = -gA2 * r_d145 /(32.0*pi^2 *Fpi4)            
            tmp_s += Lq * (k2 +3.0/8.0 *q2 + nd_mpi4 /w2) *f_N3LO_Vs
            tmp_s += w2 * Lq * f_N3LO_2l_Vs
        end
    end
    ## N4LO
    if chi_order >= 4 && EMN && calc_TPE_separately
        obj = LoopObjects.n4lo; tmp_s += Integral_ImV_T(obj.mudomain,q2,obj,obj.ImVt)  # VsB(2.15) 
        obj = LoopObjects.n4lo; tmp_s += Integral_ImV_T(obj.mudomain3,q2,obj,obj.ImVs3)# Vs12(2.26) Vt13(2.31) Vt14(2.35)
    end
    tmp_s *= -q2  
    return tmp_s
end

"""
    Ws_term(chi_order,LoopObjects,w,q2,Lq,Aq,nd_mpi,c4,Fpi2;EMN=false,useMachleidt=true,calc_TPE_separately=false)

- NNLO: EM Eq.(4.16) & Eq.(4.24) => `it_pi` term
- N3LO: EM Eq.(D.2) (``c^2_i`` term), Eq.(D.6) (`c_i/M_N`` term), Eq.(D.12) (`1/M^2_N`` term), Eq.(D.27) (2-loop term)
"""
function Ws_term(chi_order,LoopObjects,w,q2,Lq,Aq,nd_mpi,c4,Fpi2;EMN=false,useMachleidt=true,calc_TPE_separately=false)
    Fpi4 = Fpi2^2; Fpi6 = Fpi2^3; w2 = w^2; w2Aq = w2 * Aq
    nd_mpi2 = nd_mpi^2; nd_mpi4 = nd_mpi^4
    tmp_s = 0.0
    if chi_order >= 2
        f_NNLO_Ws = - gA2 / (32.0 * pi * Fpi4) 
        if EMN 
           tmp_s += (c4 * w2Aq) * f_NNLO_Ws             
            if chi_order >=3 #leading relativistic correction 
                if useMachleidt
                    tmp_s += -1/4 * (gA2*(3*nd_mpi2 + q2) -w2) *Aq *f_NNLO_Ws
                else
                    tmp_s += -1/4 * (gA2*(3*nd_mpi2 + q2) -w2) *Aq *f_NNLO_Ws
                    tmp_s += f_NNLO_Ws * gA2 * nd_mpi /8.0
                end
            end
        else
            term1 = (c4+0.25) * w2Aq
            term2 = -0.125*gA2 * (10.0 *nd_mpi2 + 3.0 *q2) *Aq
            it_pi = gA2 /8.0 * (nd_mpi + w2Aq) #nd_mpi term is different from EMN
            tmp_s += (term1+term2+it_pi) * f_NNLO_Ws
        end
    end
    if chi_order >= 3 
        f_N3LO_a_Ws = c4^2 / (96.0*pi^2 *Fpi4)
        tmp_s += Lq * w2 * f_N3LO_a_Ws
        f_N3LO_b_Ws = - c4 / (192.0*pi^2 *Fpi4)
        tmp_s += Lq * ( gA2 *(16.0*nd_mpi2 + 7.0 * q2) -w2) * f_N3LO_b_Ws
        f_N3LO_c_Ws = 1.0/ (1536.0* pi^2 * Fpi4)
        f_N3LO_2l_Ws= gA4 / (2048.0 * pi^2 * Fpi6)
        if EMN
            #2loop
            obj = LoopObjects.n3lo; tmp_s += Integral_ImV_T(obj.mudomain,q2,obj,obj.ImWt) 
        else
            ## EM terms
            tmp_s += Lq * ( 4.0*gA4 *(7.0*nd_mpi2 +17.0/4.0 *q2 + 4.0*nd_mpi4 /w2)
                            -32.0 *gA2 * (nd_mpi2 + 7.0/16.0 * q2) + w2) * f_N3LO_c_Ws
            tmp_s += f_N3LO_2l_Ws * w2Aq *( w2Aq + 2.0*nd_mpi * (1.0+2.0*gA2))
        end
    end
    if chi_order >= 4 && EMN && calc_TPE_separately
        obj = LoopObjects.n4lo; tmp_s += Integral_ImV_T(obj.mudomain,q2,obj,obj.ImWt) # EKMN(2.13)(2.16)
        obj = LoopObjects.n4lo; tmp_s += Integral_ImV_T(obj.mudomain3,q2,obj,obj.ImWs3) # 3loop EKMN(2.28) & EKMN(2.33)
    end 
    tmp_s *= -q2        
    return tmp_s
end

"""
    Vc_term(chi_order,w,tw2,q2,Lq,Aq,nd_mpi,c1,c2,c3,Fpi2,LoopObjects;EMN=false,usingMachleidt=true,ImVerbose=false,calc_TPE_separately=false)
- NNLO: EM NNLO eq.(4.13), EKMN Eq.(C1) + relativistic correction Eq.(D7)
- N3LO: 
    - f1term ``c^2_i``: EM Eq.(D.1), EKMN Eq.(D1)
    - f6term ``c_i/M_/N``: EM Eq.(D.4)

!!! note
    While Eq.(D7) in EKMN paper are given with the statement "given by [1]" refering the EM review, Eq.(D7) doesn't match neither Eq.(4.13) nor Eq.(4.13)+Eq.(4.21).
    Since the LECs in EKMN paper are determined through the Eq.(D7) and this correction term gives minor effect, we use Eq.(D7).
    If `usingMachleidt` is set `false`, one can restore original expression in EM review.
"""
function Vc_term(chi_order,w,tw2,q2,Lq,Aq,nd_mpi,c1,c2,c3,Fpi2,LoopObjects;EMN=false,usingMachleidt=true,ImVerbose=false,calc_TPE_separately=false)
    Fpi4 = Fpi2^2; Fpi6 = Fpi2^3; w2 = w^2; tw2Aq = tw2 * Aq; nd_mpi2 = nd_mpi^2
    tmp_s = 0.0
    if chi_order >= 2                   
        f_NNLO_Vc = 3.0 * gA2 / (16.0 * pi * Fpi4)    
        if EMN
            tmp_s += - ( 2.0 * nd_mpi2 * (2.0*c1 -c3) -q2*c3) * tw2Aq * f_NNLO_Vc 
            if chi_order >= 3                
                factor = f_NNLO_Vc *gA2 * 2/16.0
                term1 = nd_mpi^5 / (2 * w2)
                term2 = tw2Aq * (q2-nd_mpi2)
                if usingMachleidt
                    tmp_s += factor * (term1+term2)
                else
                    diff = -gA2 /16.0 * (nd_mpi * w2) * f_NNLO_Vc
                    tmp_s += factor * (term1+term2) + diff 
                end
            end
        else
            term1 = gA2 * nd_mpi^5 / (16.0*w^2)
            term2 = - ( 2.0 * nd_mpi2 * (2.0*c1 -c3) -q2*(c3+3.0/16.0 * gA2)) * tw2Aq
            it_pi = -gA2 /16.0 * (nd_mpi * w2 + tw2 * tw2Aq )
            tmp_s +=  (term1+term2+it_pi) * f_NNLO_Vc
        end
    end
    if chi_order >= 3 
        f_N3LO_f1_Vc = 3.0/ (16.0 * pi^2 * Fpi4)
        brak1  = ((c2*w2)/6.0 + c3*tw2 -4.0*c1*nd_mpi2)^2 + (c2^2 *w2^2)/45.0
        tmp_s += Lq * brak1 * f_N3LO_f1_Vc                    
        ## EM eq(D.4)
        f_N3LO_f6_Vc = - gA2 / (32.0 * pi^2 * Fpi4)
        brak6  = (c2-6*c3) * q2^2 + 4.0*(6*c1+c2-3*c3) *q2 *nd_mpi2
        brak6 += 6.0*(c2-2.0*c3) * nd_mpi2^2  +  24.0 *(2.0*c1+c3) * nd_mpi2^3 / w2 
        tmp_s += Lq *( brak6 * f_N3LO_f6_Vc)
        ### 2-loop corrections
        if EMN 
           obj = LoopObjects.n3lo; ImV = obj.ImVc; tmp_s += Integral_ImV(obj.mudomain,q2,obj,ImV)
        else
            f_N3LO_2l_Vc = 3.0 * gA4 /(1024.0 *pi^2 * Fpi6)
            ## EM eq.(D.9)
            Mm2cor = Lq *(2.0*nd_mpi2^4 / w^4 + 8.0*nd_mpi2^3 / w2 -q2^2 -2.0*nd_mpi2^2 ) + nd_mpi2^3 /(2.0 *w2)
            tmp_s +=  f_N3LO_f6_Vc *gA2 *Mm2cor
            ## EM eq.(D.18)
            tmp_s += f_N3LO_2l_Vc * tw2Aq * ((nd_mpi2 + 2.0 * q2)*(2.0*nd_mpi + tw2Aq) +4.0*gA2 * nd_mpi *tw2)
        end
    end
    if chi_order >= 4 && EMN && calc_TPE_separately
        obj = LoopObjects.n4lo; tmp_s += Integral_ImV(obj.mudomain,q2,obj,obj.ImVc) #Eq.(2.12) & Eq.(2.17)
    end
    return tmp_s
end

"""
    Wc_term(chi_order,LoopObjects,w,tw2,q2,k2,Lq,Aq,nd_mpi,c4,r_d12,r_d3,r_d5,Fpi2;EMN=false,useMachleidt=true,calc_TPE_separately=false)
!!! note 
    EMN Eq.(D8) doesn't match the expression in EM review, Eq.(4.14) nor Eq(4.14)+Eq.(4.22).
    Since the LECs in EKMN paper are determined through the Eqs.(D8) and these difference gives minor effect, we use EMN eq.(D8).
"""
function Wc_term(chi_order,LoopObjects,w,tw2,q2,k2,Lq,Aq,nd_mpi,c4,r_d12,r_d3,r_d5,Fpi2;EMN=false,useMachleidt=true,calc_TPE_separately=false)
    Fpi4 = Fpi2^2; Fpi6 = Fpi2^3; nd_mpi2 = nd_mpi^2; nd_mpi4 = nd_mpi^4
    w2 = w^2; tw2Aq = tw2 * Aq 
    tmp_s = 0.0
    # NLO
    if chi_order >=1 
        f_NLO_Wc = -1.0 / (384.0 * pi^2 * Fpi4)
        brak = 4.0*nd_mpi2*(5.0 * gA4 - 4.0 * gA2 -1.0)
        brak +=  q2 * (23.0*gA4 -10.0*gA2 -1.0) +48.0* gA4 *nd_mpi4 / w2
        tmp_s += Lq * brak * f_NLO_Wc
    end
    ## NNLO EM eq.
    if chi_order >= 2 
        f_NNLO_Wc = gA2 / (128.0 * pi * Fpi4)
        term1 = 3.0* gA2 * nd_mpi^5 / w^2
        tmp_s += f_NNLO_Wc  * term1
        if EMN 
            if chi_order >=3 # leading rel. correction
                if useMachleidt
                    term2 = 2* (gA2*(3.0*nd_mpi2 +2.0*q2) - 2.0*nd_mpi^2 -q2) * tw2Aq
                    tmp_s += f_NNLO_Wc  * term2
                else
                    term2 = 2* (gA2*(3.0*nd_mpi2 +2.0*q2) - 2.0*nd_mpi^2 -q2) * tw2Aq
                    term2 += gA2 * nd_mpi * w2
                    tmp_s += f_NNLO_Wc  * term2
                end
            end
        else #EM Eq(4.14)+Eq.(4.22) nd_mpi*w2 term may be missed in EKMN paper
            term2 = - (4.0*nd_mpi^2 + 2.0*q2 -gA2*(4.0*nd_mpi2 +3.0*q2)) * tw2Aq
            it_pi = gA2 * ( nd_mpi * w2 + tw2 * tw2Aq) 
            tmp_s += (term2+it_pi) * f_NNLO_Wc
        end
    end
    ## N3LO
    if chi_order >= 3 
        #a: ci/Mn term,  EM eq.(D.5), EMN eq.(2.20)
        f_N3LO_a_Wc = -c4 / (192.0 * pi^2 * Fpi4)
        tmp_s += f_N3LO_a_Wc * (gA2 *(8*nd_mpi2+5*q2) +w2) * q2 * Lq        
        if EMN 
            # 2-loop correction EKMN eq.(D4)
            obj = LoopObjects.n3lo; tmp_s += Integral_ImV(obj.mudomain,q2,obj,obj.ImWc)
        else
            # #b: 1/Mn^2 term, EM eq.(D.10)
            f_N3LO_b_Wc = -1.0 / (768.0 * pi^2 * Fpi4)
            brak1 = Lq * ( 8.0*gA2 * (11.0/4.0 *q2^2 +5.0*nd_mpi2*q2
                                    + 3.0*nd_mpi2^2 -6.0*nd_mpi2^3 /w2
                                    -k2*(8.0*nd_mpi2 + 5.0*q2))
                        + 4.0*gA4 *(k2*(20.0*nd_mpi2 +7.0*q2 -16.0*nd_mpi2^2 / w2)
                                    +16.0 * nd_mpi2^4 /w2^2 + 12.0 *nd_mpi2^3 /w2
                                    -27.0/4.0* q2^2-11.0*nd_mpi2*q2-6.0*nd_mpi2^2)
                        +(q2-4.0*k2)*w2)
            brak2 =  16.0*gA4 *nd_mpi2^3 / w2
            tmp_s += (brak1 + brak2) * f_N3LO_b_Wc
            # 2l: two-loop correction EM eq.(D.20)
            f_N3LO_2l_Wc= 1.0/(18432.0 * pi^4 * Fpi6)
            term1 = 192.0* pi^2 *Fpi2*w2*r_d3* (2.0*gA2*tw2 -3.0/5.0 *(gA2-1.0) *w2)
            brak = 384.0 * pi^2 *Fpi2 * (tw2 * r_d12 + 4.0*nd_mpi2 * r_d5)
            brak += Lq * (4.0*nd_mpi2 *(1+2.0*gA2) + q2 * (1.0+5.0*gA2))
            brak += - (q2 * (5.0+13.0*gA2)/3.0 + 8.0 * nd_mpi2 * (1.0+2.0*gA2))
            term2 = (6.0*gA2 *tw2 - (gA2 -1.0) *w2)  * brak
            tmp_s += (term1 + term2) * Lq * f_N3LO_2l_Wc
        end
    end
    if chi_order >= 4 && EMN && calc_TPE_separately
        obj = LoopObjects.n4lo; tmp_s += Integral_ImV(obj.mudomain,q2,obj,obj.ImWc) #Eq.(2.18) 
        obj = LoopObjects.n4lo; tmp_s += Integral_ImV(obj.mudomain3,q2,obj,obj.ImWc3) #3PE Eq.(2.30) 
    end
    return tmp_s
end

"""
    Vls_term(chi_order,w,tw2,q2,Lq,Aq,nd_mpi,c2,Fpi2;EMN=false)

!!! note
    In EMN interaction, relativistic corrections are counted as N3LO.
"""
function Vls_term(chi_order,w,tw2,q2,Lq,Aq,nd_mpi,c2,Fpi2;EMN=false)
    Fpi4 = Fpi2^2; nd_mpi4 = nd_mpi^4; w2 = w^2; tw2Aq = tw2*Aq
    tmp_s = 0.0
    if chi_order >= 2
        f_NNLO_Vls = 3.0 * gA4 / (32.0 * pi * Fpi4)        
        if (EMN && chi_order >=3) || !EMN
            tmp_s += tw2Aq * f_NNLO_Vls
        end
    end
    if chi_order >= 3
        f_N3LO_a_Vls = c2 *gA2 /(8.0* pi^2 * Fpi4)
        tmp_s += f_N3LO_a_Vls * w2 * Lq
        f_N3LO_b_Vls = gA4 /(4.0* pi^2 * Fpi4)
        if !EMN
            tmp_s += f_N3LO_b_Vls * Lq * (11.0/32.0 * q2 + nd_mpi4 / w2)
        end
    end
    return tmp_s
end

function Wls_term(chi_order,w,q2,Lq,Aq,nd_mpi,c4,Fpi2;EMN=false)
    Fpi4 = Fpi2^2; nd_mpi2 = nd_mpi^2; nd_mpi4 = nd_mpi^4
    w2 = w^2; w2Aq = w2 * Aq
    tmp_s = 0.0
    if chi_order >= 2 
        if !EMN || (EMN && chi_order>=3)
            f_NNLO_Wls = gA2 * (1.0- gA2) / (32.0 *pi * Fpi4)
            tmp_s += w2Aq * f_NNLO_Wls
        end
    end
    if chi_order >= 3
        f_N3LO_a_Wls = -c4 /(48.0*pi^2 * Fpi4)
        f_N3LO_b_Wls = 1.0 /(256.0*pi^2 * Fpi4) 
        tmp_s += f_N3LO_a_Wls*Lq * (gA2 *(8.0*nd_mpi2 + 5.0*q2) + w2)
        if !EMN
            tmp_s += f_N3LO_b_Wls*Lq * (16.0*gA2 *(nd_mpi2 + 3.0/8.0 *q2)
                     + 4.0/3.0 * gA4 *(4.0*nd_mpi4 /w2 -11.0/4.0 * q2 -9.0*nd_mpi2) -w2)
        end
    end
    return tmp_s
end

function Vsl_term(chi_order,Lq,Fpi2;EMN=false)
    Fpi4 = Fpi2^2
    tmp_s = 0.0
    if chi_order >=3 && !EMN
        f_N3LO = gA4 /(32.0* pi^2 * Fpi4)
        tmp_s += f_N3LO * Lq
    end
    return tmp_s
end


"""
    tpe_for_givenJT(chiEFTobj,LoopObjects,Fpi2,tmpLECs,J,pnrank,ts,ws,fff,dwn,nd_mpi,xr,pjs_para,gis_para,opfs_para,f_idx,tVs_para,lsj,tllsj_para,tdict,V12mom,tmpsum,to;calc_TPE_sep=true)

Calculating TPE contribution in a given momentum mesh point.
"""
function tpe_for_givenJT(chiEFTobj,LoopObjects,Fpi2,tmpLECs,
                 J,pnrank,fff,dwn,nd_mpi,pjs_para,gis_para,opfs_para,
                 f_idx,tVs_para,lsj,tllsj_para,tdict,tmpsum,to;calc_TPE_sep=true)
    params = chiEFTobj.params
    LamSFR_nd = dwn * params.LambdaSFR
    ts = chiEFTobj.ts; ws = chiEFTobj.ws; xr = chiEFTobj.xr
    V12mom = chiEFTobj.V12mom
    chi_order = params.chi_order
    c1 = tmpLECs["c1"]; c2 = tmpLECs["c2"]; c3 = tmpLECs["c3"]; c4 = tmpLECs["c4"]  
    r_d12 = tmpLECs["r_d12"]; r_d3 = tmpLECs["r_d3"]; r_d5 = tmpLECs["r_d5"]
    r_d145 = tmpLECs["r_d145"]; r_e14 = tmpLECs["r_e14"]; r_e17 = tmpLECs["r_e17"]
    nd_mpi2 = nd_mpi^2
    hc3 = hc^3
    n_mesh = length(xr)
    usingSFR = ifelse(LamSFR_nd!=0.0,true,false)

    @inbounds @threads for V_i= 1:n_mesh
        x = xr[V_i]; xdwn = x*dwn; xdwn2 = xdwn^2       
        ex = sqrt(1.0+xdwn2)  
        tid = threadid()        
        gis = gis_para[tid]
        tllsj = tllsj_para[tid]#; tllsj .= org_tllsj      
        target = tmpsum[tid]
        tVs = tVs_para[tid]
        f_T,f_SS,f_C,f_LS,f_SL = opfs_para[tid]
        @inbounds for V_j = 1:n_mesh
            y = xr[V_j]; ydwn = y*dwn; ydwn2=ydwn^2
            f_sq!(f_T,xdwn,ydwn);f_ls!(f_LS,xdwn,ydwn);f_sl!(f_SL,xdwn,ydwn)
            ey = sqrt(1.0+ydwn2)
            k2=0.5*(xdwn2 + ydwn2)
            ree = 1.0/sqrt(ex*ey)
            fc = fff * hc3 * freg(x,y,2) * ree
            for i in eachindex(gis); gis[i] .= 0.0; end #gis [1:7][1:9] 
            @inbounds for n in eachindex(ts)
                tpj = @view pjs_para[tid][n,:] 
                t = ts[n]; int_w = ws[n]
                q2 = xdwn2 + ydwn2 -2.0*xdwn*ydwn*t; q = sqrt(q2)
                w2 = 4.0*nd_mpi2 + q2; w = sqrt(w2); tw2 = 2.0*nd_mpi2 + q2
                Lq,Aq = calc_LqAq(w,q,nd_mpi,usingSFR,LamSFR_nd)
       
                ## Tensor term: Vt
                tmp_s = Vt_term(chi_order,LoopObjects,w,tw2,q2,k2,Lq,Aq,nd_mpi,r_d145,Fpi2;EMN=usingSFR,calc_TPE_separately=calc_TPE_sep)
                axpy!(tmp_s*int_w,tpj,target[1])
                ## Tensor term: Wt
                tmp_s = Wt_term(chi_order,LoopObjects,w,q2,Lq,Aq,nd_mpi,c4,Fpi2;EMN=usingSFR,calc_TPE_separately=calc_TPE_sep)
                axpy!(tmp_s*int_w,tpj,target[2])
                ## sigma-sigma term: Vs
                tmp_s = Vs_term(chi_order,LoopObjects,w,tw2,q2,k2,Lq,Aq,nd_mpi,r_d145,Fpi2;EMN=usingSFR,calc_TPE_separately=calc_TPE_sep)
                axpy!(tmp_s*int_w,tpj,target[3])
                ## sigma-sigma term: Ws
                tmp_s = Ws_term(chi_order,LoopObjects,w,q2,Lq,Aq,nd_mpi,c4,Fpi2;EMN=usingSFR,calc_TPE_separately=calc_TPE_sep)
                axpy!(tmp_s*int_w,tpj,target[4])
                ## Central term: Vc
                tmp_s = Vc_term(chi_order,w,tw2,q2,Lq,Aq,nd_mpi,c1,c2,c3,Fpi2,LoopObjects;EMN=usingSFR,calc_TPE_separately=calc_TPE_sep)
                axpy!(tmp_s*int_w,tpj,target[5])
                ## Central term: Wc
                tmp_s = Wc_term(chi_order,LoopObjects,w,tw2,q2,k2,Lq,Aq,nd_mpi,c4,r_d12,r_d3,r_d5,Fpi2;EMN=usingSFR,calc_TPE_separately=calc_TPE_sep)
                axpy!(tmp_s*int_w,tpj,target[6])
                ## LS term: Vls
                tmp_s = Vls_term(chi_order,w,tw2,q2,Lq,Aq,nd_mpi,c2,Fpi2;EMN=usingSFR)
                axpy!(tmp_s*int_w,tpj,target[7])
                ## LS term: Wls
                tmp_s = Wls_term(chi_order,w,q2,Lq,Aq,nd_mpi,c4,Fpi2;EMN=usingSFR)
                axpy!(tmp_s*int_w,tpj,target[8])
                ## sigma-L term: Vsl
                tmp_s = Vsl_term(chi_order,Lq,Fpi2;EMN=usingSFR)
                axpy!(tmp_s*int_w,tpj,target[9])        
                if !calc_TPE_sep            
                    n4lo_tpe_integral!(LoopObjects,q2,int_w,tpj,target)
                end
                   
            end        

            for ch =1:9
                gi = gis[ch]
                target = tmpsum[threadid()]
                axpy!(1.0,target[ch],gi)
                target[ch] .= 0.0 
            end
            #1:Vt, 2:Wt, 3:Vs, 4:Ws, 5:Vc, 6:Wc, 7:Vls, 8:Wls, 9:Vsl
            calc_IJ_V(J,pnrank,gis[1],f_T,fc,f_idx,tVs,lsj,tllsj,tdict,V12mom,V_i,V_j,to)
            calc_IJ_V(J,pnrank,gis[2],f_T,fc,f_idx,tVs,lsj,tllsj,tdict,V12mom,V_i,V_j,to;isodep=true)
            calc_IJ_V(J,pnrank,gis[3],f_SS,fc,f_idx,tVs,lsj,tllsj,tdict,V12mom,V_i,V_j,to)
            calc_IJ_V(J,pnrank,gis[4],f_SS,fc,f_idx,tVs,lsj,tllsj,tdict,V12mom,V_i,V_j,to;isodep=true)     
            calc_IJ_V(J,pnrank,gis[5],f_C,fc,f_idx,tVs,lsj,tllsj,tdict,V12mom,V_i,V_j,to)
            calc_IJ_V(J,pnrank,gis[6],f_C,fc,f_idx,tVs,lsj,tllsj,tdict,V12mom,V_i,V_j,to;isodep =true)
            calc_IJ_V(J,pnrank,gis[7],f_LS,fc,f_idx,tVs,lsj,tllsj,tdict,V12mom,V_i,V_j,to;addtype="ls")    
            calc_IJ_V(J,pnrank,gis[8],f_LS,fc,f_idx,tVs,lsj,tllsj,tdict,V12mom,V_i,V_j,to;isodep=true,addtype="ls")
            calc_IJ_V(J,pnrank,gis[9],f_SL,fc,f_idx,tVs,lsj,tllsj,tdict,V12mom,V_i,V_j,to;addtype="sl")                
        end
    end
    return nothing
end

function n4lo_tpe_integral!(LoopObjects,q2,int_w,tpj,target)
    obj = LoopObjects.n4lo
    ts = obj.ts; ws = obj.ws
    fac1 = (obj.mudomain[2]-obj.mudomain[1])/2; fac2 = (obj.mudomain[1]+obj.mudomain[2])/2
    fac1_3l = (obj.mudomain3[2]-obj.mudomain3[1])/2; fac2_3l = (obj.mudomain3[1]+obj.mudomain3[2])/2
    fac_ImV = - 2 * q2^3 / pi 
    fac_ImVT = 2 * q2^2 /pi 
    sumVt = sumWt = sumVs = sumWs = sumVc = sumWc = 0.0
    @inbounds for (ith,t) in enumerate(ts)
        # ImV (Vc/Wc)
        mu = fac1*t + fac2
        ImVc = obj.ImVc[ith]; ImWc = obj.ImWc[ith]
        sumVc += fac1 * fac_ImV * ws[ith] * ImVc / (mu^5 *(mu^2 + q2))
        sumWc += fac1 * fac_ImV * ws[ith] * ImWc / (mu^5 *(mu^2 + q2))
 
        # ImVT (Vt/Wt/Vs/Ws)
        ImVt = obj.ImVt[ith]; ImWt = obj.ImWt[ith]
        sumVt += fac1 * fac_ImVT * ws[ith] * ImVt / (mu^3 * (mu^2 + q2))
        sumVs += fac1 * fac_ImVT * ws[ith] * ImVt / (mu^3 * (mu^2 + q2))
        sumWt += fac1 * fac_ImVT * ws[ith] * ImWt / (mu^3 * (mu^2 + q2))
        sumWs += fac1 * fac_ImVT * ws[ith] * ImWt / (mu^3 * (mu^2 + q2))

        ## 3-loop terms
        mu = fac1_3l * t + fac2_3l
        ImVt3 = obj.ImVt3[ith]; ImWt3 = obj.ImWt3[ith]; ImVs3 = obj.ImVs3[ith]; ImWs3 = obj.ImWs3[ith]; ImWc3 = obj.ImWc3[ith]
        sumVt += fac1_3l * fac_ImVT *  ws[ith] * ImVt3 / (mu^3 * (mu^2 + q2))
        sumWt += fac1_3l * fac_ImVT *  ws[ith] * ImWt3 / (mu^3 * (mu^2 + q2))
        sumVs += fac1_3l * fac_ImVT *  ws[ith] * ImVs3 / (mu^3 * (mu^2 + q2))
        sumWs += fac1_3l * fac_ImVT *  ws[ith] * ImWs3 / (mu^3 * (mu^2 + q2))
        sumWc += fac1_3l * fac_ImV  *  ws[ith] * ImWc3 / (mu^5 * (mu^2 + q2))
    end
    axpy!(sumVt*int_w,tpj,target[1])
    axpy!(sumWt*int_w,tpj,target[2])
    axpy!(sumVs*int_w,tpj,target[3])
    axpy!(sumWs*int_w,tpj,target[4])
    axpy!(sumVc*int_w,tpj,target[5])
    axpy!(sumWc*int_w,tpj,target[6])  

    return nothing
end
function calc_IJ_V(J,pnrank,gi,opf,fc,f_idx,tVs,lsj,tllsj,tdict,V12mom,V_i,V_j,to;isodep=false,addtype="")
    if J==0;gi[3]=gi[5]=gi[7]=0.0;end
    IJ0 = gi[1]; IJ1 = gi[2]
    IJ2 = (J* gi[2] + gi[3]) /(J+1) 
    IJ3 = sqrt(J/(J+1)) * (gi[2]-gi[3])

    V0  = opf[1] * IJ0 + opf[2] *IJ1
    V1  = opf[3] * IJ0 + opf[4] *IJ2
    V12 = opf[5] * IJ0 + opf[6] *IJ1
    V34 = opf[4] * IJ0 + opf[3] *IJ2
    V55 = opf[7] * IJ3
    V66 = opf[8] * IJ3
    IJ4 = gi[4] 
    IJ5 = (J* gi[4] + gi[5]) /(J+1) 
    IJ6 = sqrt(J/(J+1)) * (gi[4]-gi[5])
    IJ10 = gi[6]
    IJ11 = (J * gi[6] + gi[7]) / (J+1)
    IJ12 = sqrt(J/(J+1)) * (gi[6]-gi[7])
    if addtype == "ls"        
        V0  += opf[9] * IJ4 
        V1  += opf[10] * IJ1 + opf[9] * IJ5 
        V12 += opf[10] * IJ4
        V34 += opf[10] * IJ5 +opf[9] * IJ1
        e1 = opf[11] * IJ6
        V55 += e1
        V66 += e1
    elseif addtype == "sl"
        V0  += opf[6] * IJ4 
        V1  += opf[1] * IJ4 + opf[9] * IJ5 
        V12 += opf[1] * IJ10
        V34 += opf[9] * IJ1 + opf[1] * IJ11
        e1 = opf[6] * IJ12
        V55 += e1
        V66 += e1
    end        
    transV_into_lsj(J,pnrank,tVs,V0,V1,V12,V34,V55,V66;isodep=isodep)
    @inbounds for idx = 1:f_idx
        @views tllsj[2:5] .= lsj[idx]
        V12idx = get(tdict,tllsj,-1)
        if V12idx == -1;continue;end
        tfac = tVs[idx] * fc
        V12mom[V12idx][V_i,V_j] += tfac
    end
    return nothing
end

function transV_into_lsj(J,pnrank,Vs,v1,v2,v3,v4,v5,v6;isodep=false)
    ttis=1.0;ttit=1.0
    d2j1 = 1.0/(2*J+1)
    if J==0; v2=v4=v5=v6=0.0;end
    if J%2==1 && pnrank!=2; v1 =0.0;end
    if J%2==0 && pnrank!=2; v2 =0.0;end
    if J%2!=0 && pnrank!=2; v3=v4=v5=v6=0.0;end
    v34 = -sqrt(J*(J+1)) *(v3-v4)
    v56 = sqrt(J*(J+1)) * (v5+v6)
    if isodep
        is = J%2 + 1
        it = is%2 +1
        ttis = ifelse(is==2,-3.0,1.0)
        ttit = ifelse(it==2,-3.0,1.0)
    end
    Vs[1] = v1 * ttis
    Vs[2] = v2 * ttit
    Vs[3] = d2j1 * ((J+1)*v3 + J*v4-v56) * ttis
    Vs[4] = d2j1 * ( J*v3+(J+1)*v4 +v56) * ttis
    Vs[5] = -d2j1 * (v34-(J+1)*v5+J*v6) * ttis
    Vs[6] = -d2j1 * (v34+J*v5-(J+1)*v6) * ttis
    return nothing
end

struct n3lo_2loopObj
    nmu::Int64
    mudomain::Vector{Float64}
    ts::Vector{Float64}
    ws::Vector{Float64}
    ImVc::Vector{Float64}    
    ImVt::Vector{Float64} # also for ImVs
    ImWc::Vector{Float64}
    ImWt::Vector{Float64} # also for ImWs
end
struct n4lo_23loopObj
    mudomain::Vector{Float64}
    mudomain3::Vector{Float64}
    ts::Vector{Float64}
    ws::Vector{Float64}
    ImVc::Vector{Float64}    
    ImWc::Vector{Float64}
    ImVt::Vector{Float64}    
    ImWt::Vector{Float64}
    ImWc3::Vector{Float64}
    ImVt3::Vector{Float64}    
    ImWt3::Vector{Float64}
    ImVs3::Vector{Float64}
    ImWs3::Vector{Float64}
end
struct LoopObjects
    n3lo::n3lo_2loopObj
    n4lo::n4lo_23loopObj
end

"""
```math
V_{C,S}(q) = -\\frac{2q^6}{\\pi} \\int^{\\tilde{\\Lambda}}_{nm_\\pi} d\\mu 
\\frac{\\mathrm{Im} V_{C,S}(i\\mu)}{\\mu^5 (\\mu^2+q^2)} \\\\
V_T(q) = \\frac{2q^4}{\\pi}  \\int^{\\tilde{\\Lambda}}_{nm_\\pi} d\\mu 
```
"""
function precalc_2loop_integrals(chiEFTobj,LamSFR_nd,nd_mpi,Fpi2,c1,c2,c3,c4,r_d12,r_d3,r_d5,r_d145,r_e14,r_e17)
    params = chiEFTobj.params
    ts = chiEFTobj.ts; ws = chiEFTobj.ws
    mudomain = [2*nd_mpi,LamSFR_nd]; mudomain3 = [3*nd_mpi,LamSFR_nd]    
    n3loobj = def_n3lo_2loopObj(params,nd_mpi,Fpi2,r_d12,r_d3,r_d5,r_d145,mudomain,ts,ws)
    n4loobj = def_n4lo_23loopObj(params,nd_mpi,Fpi2,c1,c2,c3,c4,r_e14,r_e17,mudomain,mudomain3,ts,ws)
    return LoopObjects(n3loobj,n4loobj)
end

function def_n4lo_23loopObj(chiEFTobj,nd_mpi,Fpi2,c1,c2,c3,c4,r_e14,r_e17,mudomain2,mudomain3,ts,ws)
    TF = chiEFTobj.chi_order>=4 && chiEFTobj.pottype=="emn500n4lo"
    nmu = ifelse(TF,length(ts),0)
    ImVc = zeros(Float64,nmu); ImWc = zeros(Float64,nmu); ImVt = zeros(Float64,nmu)
    ImWt = zeros(Float64,nmu); ImWc3 = zeros(Float64,nmu); ImVt3 = zeros(Float64,nmu)
    ImWt3 = zeros(Float64,nmu); ImVs3 = zeros(Float64,nmu); ImWs3 = zeros(Float64,nmu)
    obj = n4lo_23loopObj(mudomain2,mudomain3,ts,ws,ImVc,ImWc,ImVt,ImWt,ImWc3,ImVt3,ImWt3,ImVs3,ImWs3)
    if TF
        calc_n4lo_ImVW!(obj,nd_mpi,Fpi2,c1,c2,c3,c4,r_e14,r_e17,ts,ws)
    end
    return obj
end
function calc_n4lo_ImVW!(obj::n4lo_23loopObj,mpi,Fpi2,c1,c2,c3,c4,r_e14,r_e17,ts,ws)
    n4lo_ImVW_classAB!(obj,mpi,Fpi2,c1,c2,c3,c4,r_e14,r_e17,ts,ws)    
    n4lo_ImVW_3loop!(obj,mpi,Fpi2,c1,c2,c3,c4,ts,ws)
    return nothing
end
function n4lo_ImVW_classAB!(obj,mpi,Fpi2,c1,c2,c3,c4,r_e14,r_e17,ts,ws)    
    Fpi4 = Fpi2^2; Fpi6 = Fpi2^3
    mudomain = obj.mudomain
    fac1 = (mudomain[2]-mudomain[1])/2; fac2 = (mudomain[1]+mudomain[2])/2
    Vc = obj.ImVc; Wt = obj.ImWt; Vt = obj.ImVt; Wc = obj.ImWc
    facVcA = - mpi^5 / (4096*pi^2 * Fpi6)  
    facWtA = c4 * gA2 *mpi^5 / (4096 * pi^2 * Fpi6)
    facVtB = gA4 *mpi^5 *(c3-c4) / (4096 * pi^2 * Fpi6)
    facWtB = facVcB = gA2 * mpi^5 / (4096* pi^2 * Fpi6)    
    for ith in eachindex(ts)
        t = ts[ith]; mu = fac1*t + fac2; u = mu/mpi; u2 = u^2
        u2m4 = u2-4.0
        Bu = log( (u+sqrt(u2m4))/2)
        ## ImVcA EKMN eq.(2.12)
        term1 = gA2 * sqrt(u2m4) * (5-2*u2-2/u2) * (24*c1+c2*u2m4+6*c3*(u2-2)) * log((u+2)/(u-2))
        term1 += 8/u * ( 3*(4*c1+c3*(u2-2))*(4*gA4*u2-10*gA4+1) + c2*(6*gA4*u2-10*gA4-3)) * Bu
        term2 = sqrt(u2m4) * ( 3*(2-u2)*(4*c1+c3*(u2-2)) + c2*(7*u2-6-u2^2) + 4*gA2 /u * (2*u2-1) * (4*(6*c1-c2-3*c3)+(c2+6*c3)*u2)
                                +4*gA4 *( 32/(u+2) *(2*c1+c3) + 64/(3*u) *(6*c1+c2-3*c3) +14*c3-5*c2-92*c1 +8*u*(18*c3-5*c2)/3
                                + u2/6 * (36*c1+13*c2-156*c3) + u2^2 /6 *(2*c2+9*c3) ))        
        Vc[ith] += facVcA * (term1+term2)
        ## ImWtA EKMN eq.(2.13)
        tV = 8*gA2 *u*(5-u2) * Bu + 1/3 * u2m4^(5/2) *log((u+2)/(u-2)) + u*sqrt(u2m4) /3 * (gA2*(30*u-u^3-64)-4*u2+16)
        Wt[ith] += facWtA * tV / mu^2
        ## ImVtB EKMN eq.(2.15)
        tV = u * (sqrt(u2m4)*(u^3-30*u+64) +24*(u2-5)*Bu)
        Vt[ith] += facVtB * tV / mu^2
        ## ImWtB EKMN eq.(2.16)
        tV = (4-u2) * ( c4/3 * (sqrt(u2m4)*(2*u2-8)*Bu+4*u*(2+9*gA2)-5*u^3 /3) + 2*r_e17*(64* pi^2 *Fpi2) *(u^3-2*u) )
        Wt[ith] += facWtB * tV / mu^2
        ## ImVcB EKMN eq.(2.17)
        tV = facVcB * (u2-2) * (1/u2 -2) * ( 2*sqrt(u2m4)*(24*c1+c2*u2m4+6*c3*(u2-2))*Bu + u*(c2*(8-5*u2/3)+6*c3*(2-u2)-24*c1))
        tV += 3*gA2*(mpi^5) *((2-u2)^3) *r_e14 / (16 *Fpi4 * u) 
        Vc[ith] += tV 
        ## ImWc EKMN eq.(2.18)
        tV  = -c1 * mpi^5 / (64*Fpi6*pi^2) * ( (3*gA2+1)/8 *sqrt(u2m4) *(2-u2) + ( (3*gA2+1)/u - 2*gA2*u)*Bu)
        tV += -c2 * mpi^5 / (64*Fpi6*pi^2) * (sqrt(u2m4)/96 *(7*u2-6-u2^2+gA2*(5*u2-6-2*u2^2)) + (gA2*u2-1-gA2)*Bu/(4*u))
        brak1 = 2*sqrt(u2m4)/9 *( 3*(7*u2-6-u2^2) + 4*gA2 *(32/u-12-20*u+7*u2-u2^2) + gA4*(114-512/u+368*u-169*u2+7*u2^2+192/(u+2)))
        brak2 = 16/(3*u) * (gA4*(6*u2^2-30*u2+35)+gA2*(6*u2-8)-3)*Bu
        tV += -c3 * mpi^5 / (4096*Fpi6*pi^2) * ( brak1 + brak2)
        brak1 = 2*sqrt(u2m4)/9 * (30-128/u +80*u -13*u2 -2*u2^2 +gA2*(512/u-114-368*u+169*u2-7*u2^2-192/(u+2)))
        brak2 = 16/(3*u) * (5-3*u2+gA2*(30*u2-35-6*u2^2))*Bu
        tV += -c4 * gA2 *mpi^5 /(4096*Fpi6*pi^2) * (brak1+brak2) 
        Wc[ith] += tV        
    end
    return nothing
end
function n4lo_ImVW_3loop!(obj,mpi,Fpi2,c1,c2,c3,c4,ts,ws)    
    Fpi6 = Fpi2^3; mudomain = obj.mudomain3
    fac1 = (mudomain[2]-mudomain[1])/2; fac2 = (mudomain[1]+mudomain[2])/2
    Wc=obj.ImWc3; Vt=obj.ImVt3; Wt=obj.ImWt3; Vs=obj.ImVs3; Ws=obj.ImWs3
    for ith in eachindex(ts)
        t = ts[ith]; mu = fac1*t + fac2; u = mu/mpi; u2 = u^2; u3= u2*u; u4=u3*u; u5=u4*u; u6=u5*u; mpi5 = mpi^5
        y = sqrt((u-3)*(u+1)); Du = log( (u-1+y)/2)
        ## ImWc13 EKMN eq.(2.30)
        Wc[ith] += -gA4 * c4 *mpi5 / (4096*Fpi6*pi^2) * (8*y/3 *(u-1)*(u-4-2*u2-u3) + 32*Du *(u3-4u+1/u))
        ## ImVs/ImVt EKMN eq.(2.31)
        tV = -gA4 * c4 *mpi5 / (4096*Fpi6*pi^2 * u3) * (y/24 *(u-1)*(37*u2^3+74*u5 -251*u4-268*u3+349*u2-58*u-135)+2*Du*(39*u4-2-52*u2-6*u2^3)) / mu^2
        Vs[ith] += tV; Vt[ith] += tV  
        Vt[ith] += -gA4*c4*mpi^3 / (4096*Fpi6*pi^2 *u5) * ( y/12 *(u-1)*(5*u2^3+10*u5-3*u4-252*u3-443*u2-58*u-135)+4*Du*(3*u4+22*u2-2))
        ## ImWs/ImWt EKMN eq.(2.33)
        brak1 = 2*c1*u*(5*u3+10*u2-5*u-4) +c2/48 * (135+58*u-277*u2-36*u3+147*u4-10*u5-5*u2^3)
        brak1 += c3/8 *(7*u2^3+14*u5-145*u4-20*u3+111*u2+18*u+27)
        brak1 += c4/6 *(44*u3+37*u4-14*u5-7*u2^3-3*u2-18*u-27)
        tV = -gA4* mpi5 / (4096*Fpi6*pi^2 *u3 * mu^2)  * ( y*(u-1)*brak1 + Du*(24*c1*(1+4*u2-3*u4)+c2*(2+2*u2-3*u4) +6*c3*u2*(3*u2-2)+8*c4*u2*(u4-5*u2+5)))
        Ws[ith] += tV; Wt[ith] += tV
        brak = 4*c1*u*(5*u3+10*u2+7*u-4) + c2/24 *(135+58*u+227*u2+204*u3+27*u4-10*u5-5*u6)
        brak += c3/4 *(27+18*u-9*u2-68*u3-121*u4+14*u5+7*u6) + c4*(4*u3+19*u4-2*u5-u6-9*u2-6*u-9)
        Wt[ith] += -gA4*mpi^3 / (4096*Fpi6*pi^2 *u5) * ( y*(u-1) *brak + 2*Du* (24*c1*(1-3*u4)+c2*(2-10*u2-3*u4)+6*c3*u2*(3*u2+2)-8*c4*u4))
        ## ImVs14 EKMN eq.(2.35)
        Vs[ith] += -gA4*c4*mpi5 /(4096*Fpi6*pi^2 *u3) * (y/24 *(u-1)*(637*u2-58*u-135+116*u3-491*u4-22*u5-11*u6) +2*Du*(6*u6-9*u4+8*u2-2) )/ mu^2
        ## ImVs12 EKMN eq.(2.26) / ImVt12 EKMN eq.(2.27)
        tV= -gA2*c4*mpi5/(4096*Fpi6*pi^2 *u3 *mu^2) * (y/12*(u-1)*(100*u3-27-50*u-151*u2+185*u4-14*u5-7*u6)+4*Du*(2+10*u2-9*u4))
        Vs[ith] += tV; Vt[ith] += tV
        Vt[ith] += -gA2*c4*mpi^3 /(4096*Fpi6*pi^2 *u5) * (y/6*(u-1)*(u6+2*u5-39*u4-12*u3+65*u2-50*u-27)+8*Du*(3*u4-10*u2+2))
        ## ImWs12 EKMN eq.(2.28)
        brak = 4*c1*u/3*(u3+2*u2-u+4)+c2/72*(u6+2*u5-39*u4-12*u3+65*u2-50*u-27)
        brak += c3/12*(u6+2*u5-31*u4+4*u3+57*u2-18*u-27)+c4/72*(7*u6+14*u5-185*u4-100*u3+151*u2+50*u+27)
        tV = -gA2*mpi5 /(4096*Fpi6*pi^2 *u3 *mu^2) * ( y*(u-1)*brak + Du*(16*c1*(4*u2-1-u4)+2*c2/3*(2-10*u2+3*u4)+4*c3*u2*(u2-2)+2*c4/3*(9*u4-10*u2-2)))
        ## ImWt12 EKMN eq.(2.29)
        Ws[ith] += tV; Wt[ith] += tV
        brak = 16*c1*u/3*(2+u-2*u2-u3) + c2/36*(73*u4-6*u5-3*u6+44*u3-43*u2-50*u-27) 
        brak += c3/2 *(19*u4-2*u5-u6+4*u3-9*u2-6*u-9) +c4/36*(39*u4-2*u5-u6+12*u3-65*u2+50*u+27)
        Wt[ith] += - gA2*mpi^3 / (4096*Fpi6*pi^2 *u5) * (y*(u-1)*brak + 4*Du*(8*c1*(u4-1)+c2*(2/3-u4)-2*c3*u4+c4/3*(10*u2-2-3*u4)))
    end
    return nothing
end
   
function def_n3lo_2loopObj(chiEFTobj,nd_mpi,Fpi2,r_d12,r_d3,r_d5,r_d145,mudomain,ts,ws)
    TF = chiEFTobj.chi_order>=3 && occursin("emn",chiEFTobj.pottype)
    nmu = ifelse(TF,length(ts),0)
    ImVc = zeros(Float64,nmu); ImVt = zeros(Float64,nmu)
    ImWc = zeros(Float64,nmu); ImWt = zeros(Float64,nmu)    
    obj = n3lo_2loopObj(nmu,mudomain,ts,ws,ImVc,ImVt,ImWc,ImWt)
    if TF 
        calc_n3lo_ImVW!(obj,nd_mpi,Fpi2,r_d12,r_d3,r_d5,r_d145,mudomain,ts,ws) 
    end
    return obj
end
function calc_n3lo_ImVW!(obj::n3lo_2loopObj,mpi,Fpi2,r_d12,r_d3,r_d5,r_d145,mudomain,ts,ws)
    n3lo_ImVc!(obj.ImVc,mpi,Fpi2,mudomain,ts,ws)
    n3lo_ImWc!(obj.ImWc,mpi,Fpi2,r_d12,r_d3,r_d5,r_d145,mudomain,ts,ws)
    n3lo_ImVsVt!(obj.ImVt,mpi,Fpi2,r_d145,mudomain,ts,ws)
    n3lo_ImWsWt!(obj.ImWt,mpi,Fpi2,mudomain,ts,ws)
    return nothing
end
"""
    n3lo_ImVc!(muvec,V,mpi,Fpi)
The expression can be found in eq.(D3) of EKMN paper.
"""
function n3lo_ImVc!(V,mpi,Fpi2,mudomain,ts,ws)
    mpi2 = mpi^2; Fpi6 = Fpi2^3
    fac1 = (mudomain[2]-mudomain[1])/2; fac2 = (mudomain[1]+mudomain[2])/2    
    prefac = 3 * gA4 / (pi*4096*Fpi6)
    @threads for ith in eachindex(ts)
        t = ts[ith]; mu = fac1*t + fac2; mu2 = mu^2
        brak = (mpi2-2*mu2) * (2*mpi +  (2*mpi2 - mu2)/(2*mu) *log( (mu+2*mpi)/(mu-2*mpi))) 
        brak += 4*gA2*mpi * (2*mpi2-mu2)       
        fac = (2*mpi2 - mu2) / mu 
        V[ith] = prefac * fac * brak
    end
    return nothing
end

"""
    n3lo_ImWc!(V,mpi,Fpi,r_d12,r_d3,r_d5,r_d145,mudomain,ts,ws)
The expression can be found in eq.(D4) of EKMN paper.
"""
function n3lo_ImWc!(V,mpi,Fpi2,r_d12,r_d3,r_d5,r_d145,mudomain,ts,ws)
    mpi2 = mpi^2; Fpi6 = Fpi2^3
    fac1 = (mudomain[2]-mudomain[1])/2; fac2 = (mudomain[1]+mudomain[2])/2
    factor = 2/(3*512*pi^3 * Fpi6)    
    @threads for ith in eachindex(ts)
        t = ts[ith]; mu = fac1*t + fac2; mu2 = mu^2 
        kappa = sqrt(mu2/4 -mpi2)
        integ_x = n3lo_integ_x_Wc(ts,ws,mu,kappa,mpi,r_d12,r_d3,r_d5,r_d145,Fpi2)
        tV = factor * kappa /mu  * integ_x
        V[ith] = tV 
    end
    return nothing
end
function n3lo_integ_x_Wc(ts,ws,mu,kappa,mpi,r_d12,r_d3,r_d5,r_d145,Fpi2)
    fac_quad = 0.5; mpi2 = mpi^2; mu2 = mu^2
    integ_x = 0.0
    for ith in eachindex(ts)
        w = ws[ith]; t = ts[ith]; x = 0.5 * t + 0.5
        kx = kappa * x; kx2 = kx^2
        brak = gA2 *(mu2-2*mpi2) + 2*(1-gA2)* kx2
        term1 = 96* pi^2 * Fpi2 * ( (2*mpi2-mu2)*r_d12 - 2*kx2*r_d3 + 4*mpi2*r_d5)
        term2 = (4*mpi2*(1+2*gA2)-mu2*(1+5*gA2)) * kappa/mu * log( (mu+2*kappa)/(2*mpi) )
        term3 = mu2/12 * (5+13*gA2) -2 *mpi2*(1+2*gA2)
        term4 = -3*kx2 + 6 * sqrt(kx2 * (mpi2+kx2)) * log( (kx+sqrt(mpi2+kx2))/mpi)
        term5 = gA4 *(mu2-2*kx2-2*mpi2) *(5/6 + mpi2/kx2 -(1+mpi2/kx2)^(3/2) * log( (kx+sqrt(mpi2+kx2))/mpi))
        integ_x += fac_quad * w * brak * (term1+term2+term3+term4+term5)
    end
    return integ_x
end

"""
    n3lo_ImVsVt!(Vt,mpi,Fpi,r_d145,mudomain,ts,ws)
The expression can be found in eq.(D5) of EKMN paper and the second term in eq.(D5) is ommited.
Note that LECs ``\\bar{d}_{14}-\\bar{d}_{15}`` is given as d145 in NuclearToolkit.jl.
`Vt` can be used for `Vs=-q2*Vt`
"""
function n3lo_ImVsVt!(Vt,mpi,Fpi2,r_d145,mudomain,ts,ws)
    mpi2 = mpi^2; Fpi4 = Fpi2^2; Fpi6 = Fpi2^3
    fac1 = (mudomain[2]-mudomain[1])/2; fac2 = (mudomain[1]+mudomain[2])/2
    fac_integ_x_Vt = 2*(gA2^3)/(512*pi^3 * Fpi6)
    @threads for ith in eachindex(ts)
        t = ts[ith]; mu = fac1*t + fac2; mu2 = mu^2 
        kappa = sqrt(mu2/4 -mpi2)
        tV = gA2 *mu * kappa^3 * (-r_d145) / ( 8*pi*Fpi4 )
        tV += fac_integ_x_Vt * mu * kappa^3 * n3lo_integ_x_Vt(ts,ws,kappa,mpi)
        Vt[ith] = tV / mu2
    end
    return nothing
end
function n3lo_integ_x_Vt(ts,ws,kappa,mpi)
    fac_quad = 0.5
    kappa2 = kappa^2
    integ_x = 0.0
    for ith in eachindex(ts)
        w = ws[ith]; t = ts[ith]; x = 0.5 * t + 0.5; x2 = x^2
        mkx = mpi^2 /(kappa2*x2)
        nume = kappa * x + sqrt( mpi^2 + kappa2*x2)
        integ_x += fac_quad * w * (1.0-x2) * ( 1/6 - mkx + (1+mkx)^(3/2) * log(nume/mpi))
    end
    return integ_x
end

"""
    n3lo_ImWsWt!(Wt,mpi,Fpi,mudomain,ts,ws)
The expression can be found in eq.(D6) of EKMN paper.
`Wt` can be used for `Ws=-q2*Wt`
"""
function n3lo_ImWsWt!(Wt,mpi,Fpi2,mudomain,ts,ws)
    mpi2 = mpi^2; Fpi6 = Fpi2^3
    fac1 = (mudomain[2]-mudomain[1])/2; fac2 = (mudomain[1]+mudomain[2])/2
    @threads for ith in eachindex(ts)
        t = ts[ith]; mu = fac1*t + fac2; mu2 = mu^2 
        fac = gA4 * (4*mpi2-mu2) / ( pi* 4096*Fpi6)
        brak = (mpi2 - mu2/4) * log( (mu+2*mpi)/(mu-2*mpi) ) +(1+2*gA2)*mu*mpi
        tW = fac * brak
        Wt[ith] = tW / mu2
    end
    return nothing
end

function Integral_ImV(mudomain,q2,obj,ImV;verbose=false)    
    if length(ImV)==0;return 0.0;end
    fac1 = (mudomain[2]-mudomain[1])/2; fac2 = (mudomain[1]+mudomain[2])/2
    ts = obj.ts; ws = obj.ws
    valV = 0.0
    fac = - 2 * q2^3 / pi
    @inbounds for (ith,t) in enumerate(ts)
        mu = fac1*t + fac2
        valV += ws[ith] * ImV[ith] / (mu^5 *(mu^2 + q2))
    end
    return fac1 * valV * fac
end
function Integral_ImV_T(mudomain,q2,obj,ImV;verbose=false)
    if length(ImV)==0;return 0.0;end
    fac1 = (mudomain[2]-mudomain[1])/2; fac2 = (mudomain[1]+mudomain[2])/2
    ts = obj.ts; ws = obj.ws
    valV = 0.0
    fac = 2 * q2^2 /pi
    for (ith,t) in enumerate(ts)
        mu = fac1 * t + fac2
        valV += ws[ith] * ImV[ith] / (mu^3 * (mu^2 + q2))
    end
    return fac1 * fac * valV    
end
