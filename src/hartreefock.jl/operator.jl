"""
    aOp!(Op::Operator,a::Float64)
function to multiply scaler to an operator.
"""
function aOp!(Op::Operator,a)
    Op.zerobody[1] *= a 
    for pn =1:2
        lmul!(a,Op.onebody[pn])
    end
    for ch in eachindex(Op.twobody)
        lmul!(a,Op.twobody[ch])
    end
    return nothing
end

"""
    aOp1_p_bOp2!(Op1::Operator,Op2::Operator,a::Float64,b::Float64)
function to overwrite ```Op2``` by ```a*Op1 + b*Op2```
"""
function aOp1_p_bOp2!(Op1::Operator,Op2::Operator,a,b)
    aOp!(Op2,b)
    Op2.zerobody[1] += a * Op1.zerobody[1]
    for pn = 1:2
        axpy!(a,Op1.onebody[pn],Op2.onebody[pn])
    end
    Op2_2b = Op2.twobody
    for ch in eachindex(Op2_2b)
        axpy!(a,Op1.twobody[ch],Op2_2b[ch])
    end
    return nothing
end

"""
    aOp1_p_bOp2_Op3!(Op1::Operator,Op2::Operator,Op3::Operator,a,b,c)
function to overwrite `Op3` by `c*Op3 + a*Op1 + b*Op2`
"""
function aOp1_p_bOp2_Op3!(Op1::Operator,Op2::Operator,Op3::Operator,a,b,c)
    aOp!(Op3,c)
    Op3.zerobody[1] += a * Op1.zerobody[1] + b * Op2.zerobody[1]
    for pn = 1:2
        axpy!(a,Op1.onebody[pn],Op3.onebody[pn])
        axpy!(b,Op2.onebody[pn],Op3.onebody[pn])
    end
    Op3_2b = Op3.twobody
    for ch in eachindex(Op3_2b)
        axpy!(a,Op1.twobody[ch],Op3_2b[ch])
        axpy!(b,Op2.twobody[ch],Op3_2b[ch])
    end
    return nothing
end

"""
    InitOp(Chan1b,Chan2b)

initialize scaler Operator, note that hermite is true as default
"""
function InitOp(Chan1b,Chan2b)
    onebody = Matrix{Float64}[ ]
    for pn = 1:2
        t1b = Chan1b.chs1b[pn]
        tkeys = keys(t1b)
        dim = sum( [ length(tkey) for tkey in tkeys] ) 
        push!(onebody,zeros(Float64,dim,dim))
    end
    twobody = Matrix{Float64}[ ]
    nch = length(Chan2b)
    for ch = 1:nch
        tkets = Chan2b[ch].kets
        npq = length(tkets)
        push!(twobody,zeros(Float64,npq,npq))
    end
    return Operator([0.0],onebody,twobody,true,false)
end

"""
function to redefine RCM. Note that `non0_ij` for Calculate_RCM assumed to be `false`.
"""
function difA_RCM(Op::Operator,Aold,Anew)
    for pn = 1:2
        Op.onebody[pn] .*= Aold^2 / (Anew^2)
    end
    O2b = Op.twobody
    nch = length(O2b)     
    for ch = 1:nch
        O2b[ch] .*= Aold^2 / (Anew^2)
    end
    return nothing
end
function difA_TCM(Op::Operator,Aold,Anew)
    for pn = 1:2
        Op.onebody[pn] .*= Aold / (Anew)
    end
    return nothing
end

"""
TCM/A

One doesn't have to calc. two-body part here, since the necessary part is already calculated in Vpp.
"""
function CalculateTCM!(Op,binfo,Chan1b,Chan2b,sps)
    hw = binfo.hw
    A = binfo.nuc.A
    for pn = 1:2
        chan1b = Chan1b.chs1b[pn]
        onebody = Op.onebody[pn]
        for a in keys(chan1b)
            oa =sps[a]; na = oa.n; la =oa.l
            idx_a = div(a,2) + a%2
            for b in chan1b[a]
                ob =sps[b]; nb = ob.n; lb =ob.l
                idx_b = div(b,2) + b%2
                tij = 0.0
                if na == nb
                  tij = 0.5 * (2*na+la+3/2) *hw / A 
                elseif na == nb-1
                  tij = 0.5 * sqrt( nb * (nb+lb+1/2) ) *hw/A
                end
                onebody[idx_b,idx_a] = onebody[idx_a,idx_b] = tij/A
            end
        end
    end

    return nothing
end

"""
    Calculate_RCM(binfo,Chan1b,Chan2b,sps,Op_Rp2,d9j,HOBs,to;non0_cm=true,non0_ij=true)
calculate ``R_{CM}`` term

Note that rirj term is also included here to avoid redundancy.
"""
function Calculate_RCM(binfo,Chan1b,Chan2b,sps,Op_Rp2,d9j,HOBs,to;non0_cm=true,non0_ij=true)
    b2 = hc2 / (Mm * binfo.hw)
    A = binfo.nuc.A; Z = binfo.nuc.Z
    ## one-body part
    for pn = 1:2
        chan1b = Chan1b.chs1b[pn]
        onebody = Op_Rp2.onebody[pn]
        for a in keys(chan1b)
            oa =sps[a]; na = oa.n; la =oa.l
            idx_a = div(a,2) + a%2
            onebody[idx_a,idx_a] = 2*na + la + 3/2
            for b in chan1b[a]
                ob =sps[b]; nb = ob.n; lb =ob.l
                idx_b = div(b,2) + b%2
                if na == nb+1
                  onebody[idx_a,idx_b] = -sqrt( na * (na+la+1/2) )
                elseif  na == nb-1
                  onebody[idx_a,idx_b] = -sqrt( nb * (nb+lb+1/2) )
                end
                onebody[idx_b,idx_a] = onebody[idx_a,idx_b]
            end
        end
        onebody .*= b2 / A^2
    end
    ## two-body part 
    twobody = Op_Rp2.twobody
    frirj = - 4/(A*Z)
    nch = length(Chan2b)
    @timeit to "RCM 2body" @threads for ch in eachindex(Chan2b)
        tmp = Chan2b[ch]
        Tz = tmp.Tz; J=tmp.J;kets = tmp.kets
        factor_rirj = ifelse(Tz==0, 0.5, 1.0)
        if Tz > 0 ;factor_rirj=0.0;end
        npq = length(kets)
        Mat = twobody[ch]
        for ib = 1:npq
            bra = kets[ib]
            for ik = ib:npq
                ket = kets[ik]
                r1r2 = calc_single_r1r2(bra,ket,sps,J,d9j,HOBs,b2,to)                
                ## RCM term
                if non0_cm
                    tRCM = 2.0 * r1r2/(A^2)
                    if tRCM==0.0;continue;end
                    Mat[ib,ik] = tRCM
                    if ib!=ik; Mat[ik,ib] = tRCM; end
                end
                ## rirj term
                if non0_ij
                    tmp = r1r2 * factor_rirj *frirj
                    if tmp == 0.0; continue;end
                    Mat[ib,ik] += tmp
                    if ib!=ik; Mat[ik,ib] += tmp;end
                end
            end                
        end
    end
    return nothing 
end

"""
    Calculate_intR2p(binfo,Chan1b,HFobj,Op_Rp2)

calculate a part of squared point proton radius R2p.
Note that this function overwrites the onebody part of given `Op_Rp2`.
"""
function Calculate_intR2p(binfo,Chan1b,HFobj,Op_Rp2)
    MS = HFobj.modelspace
    sps = MS.sps
    b2 = hc2 / (Mm * binfo.hw)
    A = binfo.nuc.A; Z = binfo.nuc.Z
    fac = b2 *(A-2)/(A*Z)
    p_chan1b = Chan1b.chs1b[1]
    onebody = Op_Rp2.onebody[1]
    for a in keys(p_chan1b)
        oa = sps[a]; na = oa.n; la = oa.l
        idx_a = div(a,2) + ifelse(a%2==0,0,1)
        onebody[idx_a,idx_a] += (2*na + la + 3/2) * fac
        for b in p_chan1b[a]
            ob = sps[b]; nb = ob.n; lb = ob.l
            idx_b = div(b,2) + ifelse(b%2==0,0,1)
            if na == nb+1
                onebody[idx_a,idx_b] += -sqrt( na * (na+la+1/2) ) * fac 
            elseif  na == nb-1
                onebody[idx_a,idx_b] += -sqrt( nb * (nb+lb+1/2) ) * fac
            end
            onebody[idx_b,idx_a] = onebody[idx_a,idx_b] 
        end
    end
    return nothing
end

"""
    Calculate_SOterm(binfo,Chan1b,HFobj,Op_Rp2)

Calculate Spin-Orbit Correction term for Rp2.
We are using "simple expression for the correction to the mean-square charge radius due to the spin-orbit term"
in the reference below:
```math
\\langle r^2 \\rangle_{SO} =-\\frac{1}{Z}\\sum_i \\frac{\\mu_i}{M^2} (\\kappa_i+1),
```
where ``\\mu_p = 2.793 \\mu_N``, ``\\mu_n = −1.913\\mu_N``, and `` \\kappa = \\ell (\\mathrm{for } j=\\ell-1/2), -(\\ell+1) (\\mathrm{for} j=\\ell+1/2)``

Ref: A.Ong, J.C.Berengut, and V.V.Flambaum, [Phys. Rev. C 82, 014320 (2010)](https://doi.org/10.1103/PhysRevC.82.014320).
"""
function Calculate_SOterm(binfo,Chan1b,HFobj,Op_Rp2)
    MS = HFobj.modelspace; sps = MS.sps
    M2 = Mm^2 / hc2 
    Z = binfo.nuc.Z
    fac = -1.0/Z /M2
    mup = 2.793  #in mu_N
    mun = -1.913 #in mu_N
    for pn = 1:2
        chan1b = Chan1b.chs1b[pn]
        onebody = Op_Rp2.onebody[pn]
        mu = ifelse(pn==1,mup,mun)
        for a in keys(chan1b)
            idx_a = div(a,2) + a%2
            oa = sps[a]; la = oa.l; ja = oa.j
            kappa = ifelse(ja==2*la-1,la,-(la+1))
            onebody[idx_a,idx_a] += fac * (kappa+1) * mu
        end
    end
    return nothing
end

"""
    Calc_Expec(binfo,Chan1b,Chan2b,HFobj,Op_Rp2,dict_2b_ch,dict6j,MatOp,to;hfmbptlevel=true,verbose=false)

Calculate expectation value of Rp2 and its HFMBPT corrections.

Details about HFMBPT correction can be found in 
Many-Body Methods in Chemistry and Physics by Isaiah Shavitt and Rodney J. Bartlett (2009, Cambridge Molecular Science) 
or Appendix in [T. Miyagi et al., Phys. Rev. C 105, 0143022 (2022)](https://doi.org/10.1103/PhysRevC.105.014302).
"""
function Calc_Expec(binfo,Chan1b,Chan2b,HFobj,Op_Rp2,dict_2b_ch,dict6j,MatOp,to;hfmbptlevel=true,verbose=false)
    MS = HFobj.modelspace; sps = MS.sps
    p_sps = MS.p_sps; Cp = HFobj.Cp; e1b_p = HFobj.e1b_p; occ_p = MS.occ_p
    n_sps = MS.n_sps; Cn = HFobj.Cn; e1b_n = HFobj.e1b_n; occ_n = MS.occ_n
    holes = MS.holes; particles = MS.particles
    Gamma = HFobj.H.twobody
    e1b = Float64[ ]
    for i = 1:length(e1b_p)
        push!(e1b,e1b_p[i])
        push!(e1b,e1b_n[i])
    end
    ## From one-body contribution    
    R2p_1b = 0.0
    for pn = 1:2
        onebody = Op_Rp2.onebody[pn]
        C = ifelse(pn==1,Cp,Cn)
        tsps = ifelse(pn==1,p_sps,n_sps)
        occ = ifelse(pn==1,occ_p,occ_n)
        rhoO = C' * onebody * C        
        for i = 1:size(C)[1]
            Nocc = (tsps[i].j +1) * occ[i,i] 
            R2p_1b += Nocc * rhoO[i,i]
        end
    end
    ## From two-body contribution
    nch = length(Chan2b)
    Omats = Op_Rp2.twobody
    R2p_2b = 0.0
    R_MP1 = 0.0
    S1 = S2 = S3 = S4 = S7 = S8 = S9 = S10 = S11 = S12 = 0.0
    for ch = 1:nch
        tbc = Chan2b[ch] #J, prty, Tz 
        J = tbc.J; Tz = tbc.Tz
        kets = tbc.kets
        nkets = length(kets)
        if nkets == 0; continue;end
        D = @view MatOp[threadid()][1:nkets,1:nkets]
        D2 = @view MatOp[threadid()+nthreads()][1:nkets,1:nkets]
        Op2b = Omats[ch]       
        for ib = 1:nkets
            i,j = kets[ib]
            phase_ij = (-1)^( div(sps[i].j+sps[j].j,2) + 1 + J)
            idx_bra1 = div(i,2) + i%2
            idx_bra2 = div(j,2) + j%2
            for ik = 1:nkets
                k,l = kets[ik]
                idx_ket1 = div(k,2) + k%2
                idx_ket2 = div(l,2) + l%2
                phase_kl = (-1)^( div(sps[k].j+sps[l].j,2) + 1 + J)

                if Tz != 0
                    C1 = ifelse(Tz<=0,Cp,Cn); C2 = ifelse(Tz>=0,Cn,Cp)
                    D[ib,ik] = C1[idx_bra1,idx_ket1] * C2[idx_bra2,idx_ket2]
                    if i!=j
                        D[ib,ik] += C1[idx_bra2,idx_ket1] * C2[idx_bra1,idx_ket2] * phase_ij
                    end
                    if i==j; D[ib,ik] *= sqrt(2.0);end
                    if k==l; D[ib,ik] /= sqrt(2.0);end
                else
                    p_idx_bra = ifelse(i%2==1,idx_bra1,idx_bra2)
                    n_idx_bra = ifelse(i%2==1,idx_bra2,idx_bra1)
                    p_idx_ket = ifelse(k%2==1,idx_ket1,idx_ket2)
                    n_idx_ket = ifelse(k%2==1,idx_ket2,idx_ket1)
                    phase = 1.0
                    phase = ifelse(i%2==0,phase_ij,1.0)           
                    phase *= ifelse(k%2==0,phase_kl,1.0)
                    D[ib,ik] = Cp[p_idx_bra,p_idx_ket] * Cn[n_idx_bra,n_idx_ket] * phase
                end

            end
        end
        BLAS.gemm!('N','N',1.0,Op2b,D,0.0,D2)        
        BLAS.gemm!('T','N',1.0,D,D2,0.0,Op2b)

        ## calc expectation value
        Nocc = 2*J + 1
        S1b = ifelse(Tz<=0,Op_Rp2.onebody[1],Op_Rp2.onebody[2])
        for ib = 1:nkets
            α, α_ = kets[ib]
            na = sps[α].occ; naa = sps[α_].occ 
            #if sps[α].occ + sps[α_].occ != 2;continue;end
            if (na ==0.0) || (naa == 0.0);continue;end
            tobe = Op2b[ib,ib] * Nocc * na * naa # needed?
            R2p_2b += tobe
        end  
        if !hfmbptlevel;continue;end
        # MBPT correction
        Gam = Gamma[ch]
        for ib = 1:nkets
            a,b = kets[ib]
            idx_a = div(a,2)+ a%2
            idx_b = div(b,2)+ b%2
            #if (sps[a].occ ==1 || sps[b].occ ==1); continue; end            
            if (sps[a].occ !=0.0 || sps[b].occ !=0.0); continue; end            
            sqab = ifelse(a==b,sqrt(2.0),1.0)
            for ik = 1:nkets
                i,j = kets[ik]
                #if sps[i].occ == 0 || sps[j].occ ==0; continue; end
                if sps[i].occ == 0.0 || sps[j].occ ==0.0; continue; end                
                sqij = ifelse(i==j,sqrt(2.0),1.0)
                sqfac = sqab * sqij
                idx_i = div(i,2)+ i%2
                idx_j = div(j,2)+ j%2
                deno = (e1b[i]+e1b[j]-e1b[a]- e1b[b]) 
                tGam = Gam[ib,ik]
                nume = tGam * Op2b[ik,ib]
                tmp = 2 * Nocc * nume / deno # 2=F_1+F_2
                R_MP1 += tmp
                # S1 = S5 term, q must be b
                for ii = 1:nkets 
                    k,q= kets[ii] 
                    #if sps[k].occ + sps[q].occ != 1;continue;end
                    if (sps[k].occ ==0.0 && sps[q].occ ==0.0) || (sps[k].occ !=0.0 && sps[q].occ !=0.0);continue;end
                    if !( q==b || k==b);continue;end
                    phase = 1.0
                    if k==b 
                        phase *= (-1)^(div(sps[k].j+sps[q].j,2)+J+1)
                        q,k = kets[ii]
                    end
                    idx_k = div(k,2)+ k%2
                    nume = - 0.5* Nocc * tGam * Gam[ii,ik] * S1b[idx_a,idx_k] * phase
                    deno = (e1b[i]+e1b[j]-e1b[a]- e1b[b]) * (e1b[k]-e1b[a])
                    S1 += nume / deno *sqfac
                end
                # S2 = S6 term G_abij G_abcj Sci 
                for ii = 1:nkets
                    c,q = kets[ii]
                    #if sps[c].occ + sps[q].occ != 1; continue;end
                    if (sps[c].occ ==0.0 && sps[q].occ == 0.0) || (sps[c].occ !=0.0 && sps[q].occ !=0.0) ; continue;end
                    if !( c==j || q==j);continue;end
                    phase = 1.0
                    if sps[c].occ !=0.0 
                        if c != j;println("errS2");exit();end                        
                        phase *= (-1)^(div(sps[c].j+sps[q].j,2)+J+1) 
                        c = q; q = kets[ii][1]
                    end
                    idx_c = div(c,2) + c%2
                    nume = 0.5 * Nocc * tGam * Gam[ib,ii] * S1b[idx_c,idx_i] *phase
                    deno = (e1b[i]+e1b[j]-e1b[a]-e1b[b]) * (e1b[i] - e1b[c])
                    S2 += nume *sqfac / deno
                end
                # S3 term
                for ii = 1:nkets
                    q,c = kets[ii]
                    if sps[q].occ + sps[c].occ != 0.0; continue;end
                    if !( c==a || q==a);continue;end
                    phase = 1.0
                    if q != a && c ==a
                        phase *= (-1)^(div(sps[q].j+sps[c].j,2)+J+1) 
                        c,q = kets[ii]
                    end
                    idx_c = div(c,2) + c%2
                    nume = 0.5 * Nocc * tGam * Gam[ii,ik] * S1b[idx_b,idx_c] *phase
                    deno = (e1b[i]+e1b[j]-e1b[a]-e1b[b]) * (e1b[i]+e1b[j]-e1b[a]-e1b[c])
                    S3 += nume  *sqfac / deno
                end
                # S4 term 
                for ii = 1:nkets
                    q,k = kets[ii]
                    if sps[q].occ * sps[k].occ ==0.0; continue;end
                    if !( q==i || k==i);continue;end
                    phase = 1.0
                    if q != i && k == i
                        phase *= (-1)^(div(sps[q].j+sps[k].j,2)+J+1) 
                        k,q = kets[ii]
                    end
                    idx_k = div(k,2) + k%2
                    nume = - 0.5*  Nocc * tGam * Gam[ib,ii] * S1b[idx_j,idx_k]
                    deno =  (e1b[i]+e1b[j]-e1b[a]-e1b[b]) * (e1b[i]+e1b[k]-e1b[a]-e1b[b])
                    S4 += nume *sqfac/deno
                end
                # S7 term
                for ii = 1:nkets
                    c,d = kets[ii]
                    if sps[c].occ + sps[d].occ != 0.0;continue;end
                    nume = 0.125*Nocc * tGam * Gam[ib,ii] * Op2b[ii,ik]
                    deno =  (e1b[i]+e1b[j]-e1b[a]-e1b[b]) * (e1b[i]+e1b[j]-e1b[c]-e1b[d])
                    S7 += nume *sqfac/deno
                end                
                # S8 term
                for ii = 1:nkets
                    k,l = kets[ii]
                    if sps[k].occ * sps[l].occ ==0.0;continue;end
                    nume = 0.125* Nocc * tGam * Gam[ik,ii] * Op2b[ib,ii]
                    deno = (e1b[i]+e1b[j]-e1b[a]-e1b[b]) * (e1b[k]+e1b[l]-e1b[a]-e1b[b])
                    S8 += nume *sqfac/deno
                end
                # S10 term
                for ii = 1:nkets
                    c,d = kets[ii]
                    if sps[c].occ +sps[d].occ != 0.0;continue;end
                    nume = 0.125*Nocc * tGam * Op2b[ib,ik] * Gam[ii,ik]
                    deno = (e1b[i]+e1b[j]-e1b[a]-e1b[b]) * (e1b[i]+e1b[j]-e1b[c]-e1b[d])
                    S10 += nume/deno
                end
                # S11 term
                for ii = 1:nkets
                    k,l = kets[ii]
                    if sps[k].occ * sps[l].occ ==0.0;continue;end
                    nume = 0.125*Nocc * tGam * Op2b[ik,ii] * Gam[ib,ii]
                    deno = (e1b[i]+e1b[j]-e1b[a]-e1b[b]) * (e1b[k]+e1b[l]-e1b[a]-e1b[b])
                    S11 += nume *sqfac/deno
                end
            end 
        end  
    end
    if hfmbptlevel
        allhs = vcat(holes[1],holes[2])
        allps = vcat(particles[1],particles[2])
        nthre = nthreads()
        keychs = [ zeros(Int64,3) for i=1:nthre]
        S912s = [ zeros(Float64,nthre) for i=1:2]
        @inbounds @threads for idx_a in eachindex(allps)
            a = allps[idx_a]
            threid = threadid()
            keych = keychs[threid]
            oa = sps[a]; la = oa.l; ja = oa.j; tz_a = oa.tz
            S9tmp = S12tmp = 0.0
            for i in allhs
                oi = sps[i]; li = oi.l; ji = oi.j; tz_i = oi.tz
                for j in allhs
                    oj = sps[j]; lj = oj.l; jj = oj.j; tz_j = oj.tz
                    tz_ij = tz_i + tz_j
                    Jmin = div(abs(ja-jj),2)
                    Jmax = div(ja+jj,2)
                    for totJ = Jmin:Jmax                           
                        if tri_check(ja,jj,totJ*2)==false;continue;end
                        Jfac = (2.0*totJ+1.0)
                        tdict6j = dict6j[totJ+1]  
                        ehole = e1b[i] +e1b[j]
                        prty_ij = (-1)^(li+lj)
                        for b in allps
                            ob = sps[b]; lb = ob.l; jb = ob.j; tz_b = ob.tz
                            tz_ab = tz_a + tz_b
                            if tz_ab != tz_ij;continue;end
                            if (-1)^(la+lb) != prty_ij;continue;end
                            if tri_check(ji,jb,totJ*2)==false;continue;end
                            keych[1] = tz_a + tz_b; keych[2] = (-1)^(la+lb)
                            v1 = vPandya(a,b,i,j,ja,jb,ji,jj,totJ,dict_2b_ch,tdict6j,Gamma,keych)
                            if v1 == 0.0;continue;end
                            v1 = v1 / (ehole - e1b[a] - e1b[b])
                            for k in allhs
                                ok = sps[k]; lk = ok.l; jk = ok.j; tz_k = ok.tz
                                prty_kb = (-1)^(lk+lb)
                                for c in allps                                  
                                    oc = sps[c]; tz_c = oc.tz; lc = oc.l; jc = oc.j
                                    if tz_i + tz_c != tz_k + tz_b;continue;end
                                    if tz_k + tz_j != tz_a + tz_c;continue;end    
                                    if (-1)^(li+lc) != prty_kb;continue;end                                
                                    if tri_check(jk,jc,totJ*2)==false;continue;end                               
                                    # S9 term  O123 = Gamma/Gamma/Omats
                                    keych[1] = tz_i + tz_c; keych[2] = prty_kb # prty_ic 
                                    v2 = vPandya(i,c,k,b,ji,jc,jk,jb,totJ,dict_2b_ch,tdict6j,Gamma,keych)
                                    if v2!=0.0
                                        keych[1] = tz_k + tz_j; keych[2] = (-1)^(lk+lj)  
                                        v3 = vPandya(k,j,a,c,jk,jj,ja,jc,totJ,dict_2b_ch,tdict6j,Omats,keych)
                                        v3 = v3 / (e1b[k] + e1b[j] -e1b[a] -e1b[c])
                                        S9tmp += - Jfac * v1 * v2 * v3 
                                    end
                                    # S12 term O123 = Gamma/Omats/Gamma
                                    keych[1] = tz_i + tz_c; keych[2] = prty_kb # prty_ic 
                                    v2 = vPandya(i,c,k,b,ji,jc,jk,jb,totJ,dict_2b_ch,tdict6j,Omats,keych)
                                    if v2!=0.0
                                        keych[1] = tz_k + tz_j; keych[2] = (-1)^(lk+lj)  
                                        v3 = vPandya(k,j,a,c,jk,jj,ja,jc,totJ,dict_2b_ch,tdict6j,Gamma,keych)
                                        v3 = v3 / (e1b[k] + e1b[j] -e1b[a] -e1b[c])
                                        S12tmp += - Jfac * v1 * v2 * v3  
                                    end
                                end
                            end
                        end
                    end
                end
            end
            S912s[1][threid] += S9tmp
            S912s[2][threid] += S12tmp
        end
        S9 = sum(S912s[1])
        S12 = sum(S912s[2])
        R_MP2 = 2*S1 + 2*S2 + S3 + S4 + 2*S7 + 2*S8 + 2*S9 + S10 + S11 + S12
        if verbose 
            Ss = [S1,S2,S3,S4,S7,S8,S10,S9,S12]
            print_vec("RSs",Ss)
            println("R2 ",R2p_1b + R2p_2b, " R_MP1 $R_MP1 MP2: $R_MP2 ")
        end
        return R2p_1b + R2p_2b, R_MP1+R_MP2
    else
        return R2p_1b + R2p_2b, 0.0
    end
end

""" 
    calc_single_r1r2(bra,ket,sps,J,dict_9j_2n,HOBs,b2,to)

Calc ``<r_1 \\cdot r_2>`` for a given 2b-channel.
- `bra`: <ab| a&b: s.p.s. (n,l,j,tz)
- `ket`: |cd> c&d: s.p.s. (n,l,j,tz)
- `dictWS`: dict of WignerSymobls
- `d9j`: 9j-symbols (array of J->S and key=[la,ja,lb,jb,L])
"""
function calc_single_r1r2(bra,ket,sps,J,dict_9j_2n,HOBs,b2,to)
    r1r2 = 0.0
    dict9j = dict_9j_2n[J+1]
    oa = sps[bra[1]]; ob = sps[bra[2]]
    oc = sps[ket[1]]; od = sps[ket[2]]
    na = oa.n; la = oa.l; ja = oa.j; tza = oa.tz
    nb = ob.n; lb = ob.l; jb = ob.j; tzb = ob.tz
    nc = oc.n; lc = oc.l; jc = oc.j; tzc = oc.tz
    nd = od.n; ld = od.l; jd = od.j; tzd = od.tz
    fab = 2*na + la + 2*nb + lb
    fcd = 2*nc + lc + 2*nd + ld
    
    for Sab = 0:1
        Scd = Sab
        tdict9j_ab = dict9j[Sab+1][div(ja,2)+1][la+1][div(jb,2)+1][lb+1]            
        tdict9j_cd = dict9j[Scd+1][div(jc,2)+1][lc+1][div(jd,2)+1][ld+1]
        for Lab = abs(la-lb):la+lb
            if !tri_check(Lab,Sab,J);continue;end
            njab = tdict9j_ab[Lab+1]
            njab *= sqrt( (2*Lab+1) *(2*Sab+1) * (ja+1) * (jb+1))
            if njab == 0.0; continue;end
            Lcd = Lab
            njcd = tdict9j_cd[Lcd+1]
            njcd *= sqrt( (2*Lcd+1) *(2*Scd+1) * (jc+1) * (jd+1) )
            if njcd == 0.0; continue; end
            for N_ab = 0:div(fab,2)                
                for Lam_ab = 0:fab-2*N_ab #Lam_ab = CoM L for a,b
                    Lam_cd = Lam_ab
                    for lam_ab =(fab-2*N_ab-Lam_ab)%2:2:fab-2*N_ab-Lam_ab #lam_ab = rel L for a,b
                        if !tri_check(Lab,Lam_ab,lam_ab);continue;end
                        asymm_factor = (abs(tza+tzc) + abs(tza+tzd) * (-1)^(lam_ab+Sab)) /2                        
                        if asymm_factor == 0.0; continue;end     
                        lam_cd = lam_ab
                        n_ab = div(fab - 2*N_ab-Lam_ab -lam_ab,2) # determined by energy conservation                        
                        ##To get HOB(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab)            
                        mosh_ab = 0.0
                        if (2*N_ab+Lam_ab > 2*n_ab+lam_ab) && (2*na+la > 2*nb+lb)
                            phase = (-1)^(Lam_ab+lb)
                            nkey1 = get_nkey_from_key6j(n_ab,N_ab,lam_ab,Lam_ab,0)
                            nkey2 = get_nkey_from_key6j(Lab,nb,na,lb,0)
                            target = HOBs[nkey1][nkey2]
                            mosh_ab = target * phase
                        else
                            nkey1 = get_nkey_from_key6j(N_ab,n_ab,Lam_ab,lam_ab,0)
                            nkey2 = get_nkey_from_key6j(Lab,na,nb,la,0)
                            mosh_ab = HOBs[nkey1][nkey2]
                        end                                                                     
                        if mosh_ab == 0.0;continue;end
                        for N_cd = max(0,N_ab-1):N_ab+1 
                            n_cd = div(fcd-2*N_cd-Lam_cd-lam_cd,2)
                            if n_cd < 0 || (n_ab!=n_cd && N_ab !=N_cd); continue;end
                            mosh_cd = 0.0
                            if (2*N_cd+Lam_cd > 2*n_cd+lam_cd) && (2*nc+lc > 2*nd+ld)
                                phase = (-1)^(Lam_cd+ld)
                                nkey1 = get_nkey_from_key6j(n_cd,N_cd,lam_cd,Lam_cd,0)
                                nkey2 = get_nkey_from_key6j(Lcd,nd,nc,ld,0)
                                target = HOBs[nkey1][nkey2]    
                                mosh_cd = target * phase
                            else
                                nkey1 = get_nkey_from_key6j(N_cd,n_cd,Lam_cd,lam_cd,0)
                                nkey2 = get_nkey_from_key6j(Lcd,nc,nd,lc,0)
                                mosh_cd = HOBs[nkey1][nkey2]
                            end   
                            if mosh_cd == 0.0; continue;end
                            r2cm = 0.0;r2rel = 0.0
                            if n_ab == n_cd 
                                if N_ab == N_cd
            	   	                r2cm = 2*N_ab + Lam_ab+1.5
                    	        elseif N_ab == N_cd + 1
                        		    r2cm = -sqrt(N_ab * (N_ab+Lam_ab+0.5))
                                elseif N_ab == N_cd -1 
                                    r2cm = -sqrt(N_cd * (N_cd+Lam_ab+0.5))
                                end
                            end
                            if N_ab == N_cd 
                                if n_ab == n_cd 
	   	                            r2rel = 2*n_ab + lam_ab+1.5
    		                    elseif n_ab == n_cd + 1
	    	                        r2rel = -sqrt(n_ab * (n_ab+lam_ab+0.5))
                                elseif n_ab == n_cd -1 
                                    r2rel = -sqrt(n_cd * (n_cd+lam_ab+0.5))
                                end              
                            end
                            prefactor = asymm_factor * njab * njcd * mosh_ab * mosh_cd
                            r1r2 += (r2cm-r2rel) * prefactor
                        end
                    end 
                end
            end
        end
    end
    ## note that delta function below including "tz"
    deno = sqrt( (1.0+delta(bra[1],bra[2])) * (1.0+delta(ket[1],ket[2])) )
    if tza+tzb == 0; deno = 1.0;end
    r1r2 *= b2 * 0.5 / deno
    return r1r2 
end 

""" 
    Calculate_Rp(binfo,Chan1b,Chan2b,HFobj,Op_Rp2,dict_9j_2n,HOBs,dict_2b_ch,dict6j,to;hfmbptlevel=true) 

To calculate squared point proton radius and its MBPT correction.
The squared point proton radius is related to the charge radius as follows
```math
R^2_{ch} = R^2_p + \\langle r^2_p \\rangle + \\frac{N}{Z} \\langle r^2_n \\rangle + \\frac{3}{4m^2_p c^4} 
+ \\langle r^2 \\rangle_{SO},
```
where ``\\langle r^2_p \\rangle = 0.769 \\mathrm{fm}^2``, ``\\langle r^2_n \\rangle = -0.116 \\mathrm{fm}^2``,
`` \\frac{3}{4m^2_p c^4} =0.033\\mathrm{fm}^2`` is the so-called Darwin-Foldy term, and the last term is Spin-Orbit correction term.
"""
function Calculate_Rp(binfo,Chan1b,Chan2b,HFobj,Op_Rp2,dict_9j_2n,HOBs,dict_2b_ch,dict6j,MatOp,to;hfmbptlevel=true)   
    @timeit to "RCM" Calculate_RCM(binfo,Chan1b,Chan2b,HFobj.modelspace.sps,Op_Rp2,dict_9j_2n,HOBs,to)
    Calculate_intR2p(binfo,Chan1b,HFobj,Op_Rp2)
    Calculate_SOterm(binfo,Chan1b,HFobj,Op_Rp2)
    @timeit to "expec" Rp,Rp_PT = Calc_Expec(binfo,Chan1b,Chan2b,HFobj,Op_Rp2,dict_2b_ch,dict6j,MatOp,to;hfmbptlevel=hfmbptlevel)
    return Rp,Rp_PT
end

"""
    eval_rch_hfmbpt(binfo,Chan1b,Chan2bD,HFobj,Op_Rp2,dict_9j_2n,HOBs,dict6j,to)

evaluate charge radii with HFMBPT
"""
function eval_rch_hfmbpt(binfo,Chan1b,Chan2bD,HFobj,Op_Rp2,dict_9j_2n,HOBs,dict6j,MatOp,to;io=stdout)
    Chan2b = Chan2bD.Chan2b; dict_2b_ch = Chan2bD.dict_ch_JPT
    tnuc = binfo.nuc
    N = tnuc.N
    Z = tnuc.Z
    DF = 0.033 
    Rp2 = 0.8775^2
    Rn2 = -0.1149
    Rpp,Rp_MP = Calculate_Rp(binfo,Chan1b,Chan2b,HFobj,Op_Rp2,dict_9j_2n,HOBs,dict_2b_ch,dict6j,MatOp,to)
    Rch2 = Rpp + Rp2 + N/Z *Rn2 + DF
    Rch = NaN 
    try 
        Rch = @sprintf("%12.6f", sqrt(Rch2))
    catch
        Rch = "NaN"
    end 
    ctxt = ""
    if Rp_MP !=0.0        
        RPT = NaN
        try 
            RPT = @sprintf("%12.6f", sqrt(Rch2+Rp_MP))
        catch 
            RPT = "NaN"
        end
        ctxt = " HF+PT => " * RPT
    end
    println(io,"   HF point proton radius ",@sprintf("%12.6f", sqrt(Rpp))," charge radius ",Rch," ",ctxt)
    getNormalOrderedO(binfo,HFobj,Op_Rp2,Chan1b,Chan2bD,dict6j,to) 
    return nothing
end

"""
    eval_rch_imsrg(binfo,Chan1b,Chan2bD,HFobj,IMSRGobj,PandyaObj,dict_9j_2n,HOBs,dictMono,dict6j,to)

evaluate charge radii with IMSRG.
"""
function eval_rch_imsrg(binfo,Chan1b,Chan2bD,HFobj,IMSRGobj,PandyaObj,dict_9j_2n,HOBs,dictMono,dict6j,MatOp,to)
    Chan2b = Chan2bD.Chan2b; dict_2b_ch = Chan2bD.dict_ch_JPT
    tnuc = binfo.nuc
    N = tnuc.N
    Z = tnuc.Z
    DF = 0.033; Rp2 = 0.8775^2; Rn2 = -0.1149
    Op_Rp2 = InitOp(Chan1b,Chan2b)
    ## HF level
    Rpp,Rp_PT = Calculate_Rp(binfo,Chan1b,Chan2b,HFobj,Op_Rp2,dict_9j_2n,HOBs,dict_2b_ch,dict6j,MatOp,to)
    getNormalOrderedO(binfo,HFobj,Op_Rp2,Chan1b,Chan2bD,dict6j,to;firstNO=true) 
    Rpp_HF = Op_Rp2.zerobody[1]
    Rch2_HF = Rpp_HF + Rp2 + N/Z *Rn2 + DF
    Rch_HF = sqrt(Rch2_HF)
    println("   HF point proton radius ",@sprintf("%12.6f", sqrt(Rpp_HF)),
            " charge radius ",@sprintf("%12.6f", Rch_HF), " => HF+PT ",@sprintf("%12.6f", sqrt(Rch2_HF+Rp_PT)))

    ##IMSRG(2) level
    ncomm = [0]
    norms = zeros(Float64,4)
    Omega = IMSRGobj.Omega
    tOmega = deepcopy(Omega); aOp!(tOmega,0.0)
    tmpOp  = deepcopy(Op_Rp2);aOp!(tmpOp,0.0)
    tmpOp2 = deepcopy(Op_Rp2);aOp!(tmpOp2,0.0)
    Nested = deepcopy(Op_Rp2);aOp!(Nested,0.0)
    nwritten = IMSRGobj.n_written_omega[1]
    for i = 1:nwritten
        read_omega_bin!(binfo,i,tOmega)
        BCH_Transform(tOmega,Op_Rp2,tmpOp,tmpOp2,Nested,ncomm,norms,Chan1b,Chan2bD,HFobj,dictMono,dict6j,PandyaObj,to)
        aOp1_p_bOp2!(tmpOp,Op_Rp2,1.0,0.0)
    end

    Rpp = tmpOp.zerobody[1]
    Rch2 = Rpp + Rp2 + N/Z *Rn2 + DF
    Rch = sqrt(Rch2)
    println("IMSRG point proton radius ",@sprintf("%12.6f", sqrt(Rpp))," charge radius ",@sprintf("%12.6f", Rch))
    return tmpOp
end

"""
    getNormalOrderedO(binfo,HFobj,targetOp,Chan1b,Chan2bD,dict6j,to;verbose=false,undo=false,OpeqH=false,firstNO=false)

NormalOrdering for a target Operator. For now, it only supports scaler operators.
"""
function getNormalOrderedO(binfo,HFobj,targetOp,Chan1b,Chan2bD,dict6j,to;verbose=false,undo=false,OpeqH=false,firstNO=false)
    Chan2b = Chan2bD.Chan2b
    dict_idx_from_chket = Chan2bD.dict_idx_from_chket
    sps = HFobj.modelspace.sps
    nch = length(Chan2b)
    Cp = HFobj.Cp; Cn = HFobj.Cn
    holes = Int64[ ]
    for (ith,o) in enumerate(HFobj.modelspace.sps)
        if o.occ != 0.0
            push!(holes,ith)
        end
    end 
    # NO0b 
    if firstNO  ## for Operator normal ordering,... here is the first time to make HF basis...?
        onebody = [ Cp' * (targetOp.onebody[1] * Cp), Cn' * (targetOp.onebody[2] * Cn)]    
        targetOp.onebody[1] .= onebody[1]
        targetOp.onebody[2] .= onebody[2]
    end
    ## 0-body from One-body 
    for (i,oi) in enumerate(sps)
        ni = oi.occ
        if ni != 0.0
            pn = 1 + div(1+oi.tz,2)
            idx_i = div(i,2) + i%2            
            targetOp.zerobody[1] += ni*(oi.j+1) * ifelse(undo,-1.0,1.0)  * targetOp.onebody[pn][idx_i,idx_i]
        end
    end
   
    ## 0-body & 1-body from 2-body
    tkey = zeros(Int64,2)
    for ch = 1:nch
        O2b = targetOp.twobody[ch]
        tkets = Chan2b[ch].kets
        J = Chan2b[ch].J
        Tz = Chan2b[ch].Tz
        hatfactor = 2*J+1
        nket = length(tkets)
        # => NO0B
        tsum = 0.0        
        for ib = 1:nket
            p,q = tkets[ib]
            if sps[p].occ * sps[q].occ ==0.0;continue;end
            tsum += sps[p].occ * sps[q].occ * O2b[ib,ib]
        end
        targetOp.zerobody[1] += tsum * hatfactor

        # => NO1B        
        tdict = dict_idx_from_chket[ch]
        for h in holes
            jh = sps[h].j; nh = sps[h].occ
            for (a,oa) in enumerate(sps)
                ja = oa.j
                if !tri_check(jh,ja,2*J);continue;end
                pn = 1 + div(1+oa.tz,2)
                idx_a = div(a,2) + a%2
                O1b = targetOp.onebody[pn] 
                phase_ah = 1.0
                if a > h
                    tkey[1] = h; tkey[2] = a
                    phase_ah *= (-1)^(div(ja+jh,2)+J+1)
                else
                    tkey[1] = a; tkey[2] = h
                end
                sqfac_ah = ifelse(a==h,sqrt(2.0),1.0)
                idx_ah = get(tdict,tkey,0)
                if idx_ah ==0; continue;end
                for b in Chan1b.chs1b[pn][a]
                    ob = sps[b]
                    jb = ob.j
                    if !tri_check(jh,jb,2*J);continue;end
                    idx_b = div(b,2) + b%2
                    phase_bh = 1.0
                    if b > h
                        tkey[1] = h; tkey[2] = b
                        phase_bh *= (-1)^(div(jb+jh,2)+J+1)
                    else
                        tkey[1] = b; tkey[2] = h
                    end
                    idx_bh = get(tdict,tkey,0)
                    if idx_bh ==0; continue;end
                    sqfac_bh = ifelse(b==h,sqrt(2.0),1.0)
                    sqfac = sqfac_ah*sqfac_bh
                    O1b[idx_a,idx_b] += nh * hatfactor * ifelse(undo,-1.0,1.0) /(ja+1.0) * O2b[idx_ah,idx_bh] *  phase_ah * phase_bh * sqfac                    
                    if idx_a != idx_b 
                        O1b[idx_b,idx_a] = O1b[idx_a,idx_b]
                    end
                end
            end
        end
    end
    return nothing
end

"""
    show_Hamil_norm(Op::Operator;tol=1.e-6,normtype="fro")

Function to show 1b/2b norm of a given Operator. It may be usuful for debug.
"""
function show_Hamil_norm(Op::Operator;tol=1.e-6,normtype="fro")
    if normtype == "fro"
        println("norm1b p ", @sprintf("%12.5e", norm(Op.onebody[1],2)))
        println("norm1b n ", @sprintf("%12.5e", norm(Op.onebody[2],2)))
        O2b = Op.twobody
        for ch = 1:length(O2b)
            tnorm = norm(O2b[ch],2)
            if tnorm > tol
                println("ch ",@sprintf("%3i",ch), @sprintf("%12.5e",tnorm))
            end
        end
    elseif normtype == "sum"
        println("norm1b p ", @sprintf("%12.5e", sum(Op.onebody[1])))
        println("norm1b n ", @sprintf("%12.5e", sum(Op.onebody[2])))
        O2b = Op.twobody
        for ch = 1:length(O2b)
            tnorm = sum(O2b[ch])
            if tnorm > tol
                println("ch ",@sprintf("%3i",ch), @sprintf("%12.5e",tnorm))
            end
        end
    end
    return nothing
end
