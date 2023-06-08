"""
    SRG(xr,wr,V12mom,dict_pwch,to)

Similarity Renormalization Group (SRG) transformation of NN interaction in CM-rel momentum space.
"""
function SRG(chiEFTobj,to) 
    if !chiEFTobj.params.srg; return nothing;end
    n_mesh = chiEFTobj.params.n_mesh
    srg_lambda = chiEFTobj.params.srg_lambda
    xr_fm = chiEFTobj.xr_fm; wr = chiEFTobj.wr
    V12mom = chiEFTobj.V12mom; dict_pwch = chiEFTobj.dict_pwch
    l1s = [0,1,2,3,4,5,6,1,1,2,3,4,5,6,0,1,2,3,4,5]
    l2s = [0,1,2,3,4,5,6,1,1,2,3,4,5,6,2,3,4,5,6,7]
    Ss  = [0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1]
    Js  = [0,1,2,3,4,5,6,0,1,2,3,4,5,6,1,2,3,4,5,6]

    ndim = n_mesh*2
    n_ode = div(ndim*(ndim+1),2)
    nthre = 1 #nthreads()
    Vs = [ zeros(Float64,ndim,ndim) for i = 1:nthre]
    Ts = [zeros(Float64,ndim,ndim) for i = 1:nthre]
    Hs = [zeros(Float64,ndim,ndim) for i = 1:nthre]
    Hts = [zeros(Float64,ndim,ndim) for i = 1:nthre]
    etas = [zeros(Float64,ndim,ndim) for i = 1:nthre]
    Rs = [zeros(Float64,ndim,ndim) for i = 1:nthre]
    tkeys = [ zeros(Int64,5) for i=1:nthre]
    sSRG = (1.0/srg_lambda)^4
    ds = 1.e-4
    srange = 0.0:ds:sSRG
    numit = length(srange)
    nch = length(l1s)
    for nth =1:nch*3
        threid = threadid()
        tkey = tkeys[threid]
        V = Vs[threid]; T = Ts[threid]; H = Hs[threid]
        Ht= Hts[threid];eta = etas[threid]; R = Rs[threid]
        iz = -1
        if nth <= nch
            nothing
        elseif nth <= 2*nch
            iz = 0
        else
            iz = 1
        end
        pnrank = iz + 2
        tdict = dict_pwch[pnrank]
        rmass = Mp * Mn / (Mp+Mn) 
        if iz == -1; rmass = 0.5 * Mp;end
        if iz ==  1; rmass = 0.5 * Mn;end
        face =  (2.0*rmass)^2 / hc4

        ich = ifelse(nth%nch!=0,nth%nch,nch)
        l1 = l1s[ich];l2 = l2s[ich];S = Ss[ich];J = Js[ich]
        if abs(iz % 2) == 1 && (l1+S+abs(iz)) % 2 == 0; continue;end
        tkey[1] = 2*iz; tkey[2]=l1; tkey[3]=l2; tkey[4]=S;tkey[5]=J
        if l1 == l2
            V12idx = tdict[tkey]; tv = V12mom[V12idx]
            sV  = @view  V[1:n_mesh,1:n_mesh]
            sT  = @view  T[1:n_mesh,1:n_mesh]
            sH  = @view  H[1:n_mesh,1:n_mesh]
            sHt = @view Ht[1:n_mesh,1:n_mesh]
            seta = @view eta[1:n_mesh,1:n_mesh]
            tR = @view R[1:n_mesh,1:n_mesh]                
            for (i,x) in enumerate(xr_fm)
                wri = wr[i]
                for (j,y) in enumerate(xr_fm)
                    sH[i,j] = sV[i,j] = tv[i,j] * x * y * sqrt(wri*wr[j])
                end
                sT[i,i] = (x *hc)^2 / (2*rmass)
                sH[i,i] += (x *hc)^2 / (2*rmass)
            end
            if norm(sH-sH',Inf) > 1.e-6; println("Tz ",@sprintf("%3i",iz)," L $l1 L' $l2 S $S J $J norm(sH-sH') ", norm(sH-sH',Inf));end
            srg_RK4(sH,sT,sHt,sV,seta,tR,sSRG,face,ds,numit,to)
            ## overwrite V12
            for (i,x) in enumerate(xr_fm)
                for (j,y) in enumerate(xr_fm) 
                    tv[i,j] = (sV[i,j]-sT[i,j]) / ( x * y * sqrt( wr[i] * wr[j]))
                end
            end
        elseif l1+2==l2
            ## l1-l1, l2-l2, l1-l2, l2-l1 part
            tkey[2]=l1; tkey[3]=l1; V12idx = tdict[tkey]; tv1 = V12mom[V12idx]
            tkey[2]=l2; tkey[3]=l2; V12idx = tdict[tkey]; tv2 = V12mom[V12idx]
            tkey[2]=l1; tkey[3]=l2; V12idx = tdict[tkey]; tv3 = V12mom[V12idx]
            tkey[2]=l2; tkey[3]=l1; V12idx = tdict[tkey]; tv4 = V12mom[V12idx]
            for (i,x) in enumerate(xr_fm)
                for (j,y) in enumerate(xr_fm) 
                    tfac = x * y * sqrt(wr[i] * wr[j])
                    V[i,j] = tv1[i,j] * tfac 
                    V[n_mesh+i,n_mesh+j] = tv2[i,j] * tfac 
                    V[i,n_mesh+j] = tv3[i,j] * tfac 
                    V[n_mesh+i,j] = tv4[i,j] * tfac 
                end
                T[i,i] = T[n_mesh+i,n_mesh+i] = (x *hc)^2 / (2*rmass)
            end
            H .= V;H .+= T
            srg_RK4(H,T,Ht,V,eta,R,sSRG,face,ds,numit,to)
            H .= V; H .-= T # Veff = H(s) - T # H is reused as Veff            
            for (i,x) in enumerate(xr_fm)
                for (j,y) in enumerate(xr_fm)
                    deno =  x * y * sqrt(wr[i]*wr[j])
                    tv1[i,j] = H[i,j] / deno 
                    tv2[i,j] = H[n_mesh+i,n_mesh+j] /  deno 
                    tv3[i,j] = H[i,n_mesh+j] / deno 
                    tv4[i,j] = H[n_mesh+i,j] / deno 
                end
            end
        end
    end
    return nothing
end

"""
    commutator(A,B,R,fac)

wrapper function to overwrite ```R``` by ```fac*(AB-BA)```
"""
function commutator(A,B,R,fac)
    BLAS.gemm!('N', 'N',  fac, A, B, 0.0, R)
    BLAS.gemm!('N', 'N', -fac, B, A, 1.0, R)
    return nothing
end

"""
    RKstep(T,Ho,eta,R,faceta,fRK,Ht)

wrapper function to calc. a Runge-Kutta (RK) step
"""
function RKstep(T,Ho,eta,R,faceta,fRK,Ht)
    BLAS.gemm!('N', 'N',  faceta, T, Ho, 0.0, eta)
    BLAS.gemm!('N', 'N', -faceta, Ho, T, 1.0, eta) # =>eta
    BLAS.gemm!('N', 'N',  1.0, eta, Ho, 0.0, R)
    BLAS.gemm!('N', 'N', -1.0, Ho, eta, 1.0, R)
    axpy!(fRK,R,Ht)
    return nothing
end
function RKstep_mul(T,Ho,eta,R,faceta,fRK,Ht)
    mul!(eta,T,Ho,faceta,0.0)
    mul!(eta,Ho,T,-faceta,1.0)
    mul!(R,eta,Ho,1.0,0.0)
    mul!(R,Ho,eta,-1.0,1.0)
    axpy!(fRK,R,Ht)    
    return nothing
end

"""
    srg_RK4(Ho,T,Ht,Hs,eta,R,sSRG,face,ds,numit,to; r_err=1.e-8,a_err=1.e-8,tol=1.e-6)

to carry out SRG transformation with RK4
"""
function srg_RK4(Ho,T,Ht,Hs,eta,R,sSRG,face,ds,numit,to;r_err=1.e-8,a_err=1.e-8,tol=1.e-6)
    Hs .= Ho
    Ht .= Hs
    for it = 1:numit
        if it ==numit; ds = sSRG -(numit-1)*ds;end
        Ho .= Hs        
        RKstep(T,Ho,eta,R,face,ds/6.0,Ht)

        Ho .= Hs; axpy!(0.5*ds,R,Ho)
        RKstep(T,Ho,eta,R,face,ds/3.0,Ht)

        Ho .= Hs; axpy!(0.5*ds,R,Ho)
        RKstep(T,Ho,eta,R,face,ds/3.0,Ht)

        Ho .= Hs; axpy!(ds,R,Ho)
        RKstep(T,Ho,eta,R,face,ds/6.0,Ht)

        Hs .= Ht
    end
    return nothing
end

