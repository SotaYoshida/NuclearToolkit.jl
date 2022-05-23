struct op_M1
    bi::Int64
    i::Int64
    f::Int64
    phase::Bool
    fac::Float64
end
struct op_E2
    bi::Int64
    i::Int64
    f::Int64
    phase::Bool
    fac::Float64
end

struct op_M1_nd 
    bl::Int64
    Nl::Int64
    br::Int64
    Nr::Int64
    phase::Bool
    fac::Float64
end
struct op_E2_nd
    bl::Int64
    Nl::Int64
    br::Int64
    Nr::Int64
    phase::Bool
    fac::Float64
end

function tri_relation(ja,jb,jc)
    TF= true
    if ja+jb < jc ; TF=false;end
    if jb+jc < ja ; TF=false;end
    if jc+ja < jb ; TF=false;end
    if abs(ja-jb) > jc ; TF=false;end
    if abs(jb-jc) > ja ; TF=false;end
    if abs(jc-ja) > jb ; TF=false;end
    return TF
end
    
function init_ho_by_mass(A,ihwprm=1)
    nmass = [938.27231, 939.56563]
    amass = (nmass[1]+nmass[2])/2
    hc = 197.32696
    if ihwprm == 1 ## 41A^(-1/3) MeV
        hw = 41.0 * (A^(-1.0/3.0))
    elseif ihwprm == 2 ## 45A^(-1/3)-25A^(-2/3) MeV
        #J. Blomqvist and A. Molinari, Nucl. Phys. A106, 545 (1968).
        hw = 45.0 * (A^(-1.0/3.0)) -25.0 * (A^(-2.0/3.0))
    end
    bpar = hc/sqrt(hw*amass)
    return hw,bpar
end

function radius_power(k,n1,l1,n2,l2,bpar)    
    s=0.0
    ll  = l1+l2+k
    ll1 = l2-l1+k
    ll2 = l1-l2+k
    ll1 = div(ll1,2)
    ll2 = div(ll2,2)

    imin = maximum([0, n1-ll1, n2-ll2])
    imax = minimum([n1, n2])
    for i = imin:imax
        s +=  binomial(ll1, n1-i) * binomial(ll2, n2-i) * (
            doublefactorial(ll+2*i+1)/doublefactorial(2*i))
    end        
    s *= sqrt(doublefactorial(2*n1)/doublefactorial(2*n1+2*l1+1)) 
    s *= sqrt(doublefactorial(2*n2)/doublefactorial(2*n2+2*l2+1))
    s = Float64(s) * (1.0/sqrt(2^k)) * (-1)^(n1-n2) * bpar^k
    return s 
end

function prep_obtr_diag(p_sps,n_sps,
                        mstates_p,mstates_n,tdims,
                        jocc_p,jocc_n,
                        pbits,nbits,bpar,
                        gfactors,eff_charge)
    ## Note that mu is in units of muN (= 0.105 efm)
    lp=length(p_sps); ln=length(n_sps); lblock = length(nbits)
    ec_p,ec_n = eff_charge
    c_mom_mu = sqrt(3.0/(4.0*pi))
    c_mom_Q = sqrt(5.0/(16.0*pi))

    ## M1 (mu; magnetic dipole moment)
    pjumps_mu = op_M1[  ] 
    njumps_mu = op_M1[  ] 
    ## E2 (Q; electric quadropole moment)
    pjumps_Q = op_E2[ ] 
    njumps_Q = op_E2[ ]

    vec_p_ani = [ false for i=1:length(mstates_p)]
    vec_p_cre = [ false for i=1:length(mstates_p)]
    ex_r2 = 0.0 ; coeff= 0.0
    for i=1:lp
        ni,li,ji,tzi = p_sps[i]
        for f = i:lp
            nf,lf,jf,tzf = p_sps[f]
            ex_r2 = radius_power(2,ni,li,nf,lf,bpar)
            l6j = wigner6j(Float64,lf,li,1, ji//2,jf//2,1//2)
            s6j = wigner6j(Float64,1//2,1//2,1,ji//2,jf//2,lf)
            Q3j = wigner3j(Float64,ji//2,2,jf//2,1//2,0,-1//2)
            c,mc_s,mc_idxs = possible_mz(p_sps[i],mstates_p)
            a,ma_s,ma_idxs = possible_mz(p_sps[f],mstates_p)

            for mode in ["M1", "E2"]
                if mode == "M1" # selection rule
                    if (ni !=nf || li != lf); continue;end
                end
                @inbounds for (ic,mc) in enumerate(mc_s)
                    vec_p_ani .= false
                    vec_p_ani[mc_idxs[ic]] = true
                    bit_c = bitarr_to_int(vec_p_ani)
                    @inbounds for (ia,ma) in enumerate(ma_s)
                        if ma != mc; continue;end #M conservation
                        vec_p_cre .= false
                        vec_p_cre[ma_idxs[ia]] = true
                        bit_a = bitarr_to_int(vec_p_cre) 
                        pfac = sqrt((ji+1)*(jf+1))
                        if mode =="M1"
                            lam = 1
                            sign1 = (-1)^((jf-ma)/2+lf+(ji+3)/2)
                            sign2 = (-1)^((jf-ma)/2+lf+(jf+3)/2)
                            fac1 = sign1 * l6j .* sqrt(lf*(lf+1)*(2*lf+1))
                            fac2 = sign2 * s6j .* sqrt(1.5)
                            m3j = pfac * wigner3j(jf//2,lam,ji//2,-ma//2,0,mc//2) 
                            coeff = m3j * (fac1*gfactors[1]+fac2*gfactors[3]) *c_mom_mu
                        elseif mode =="E2"
                            lam = 2
                            signQ = (-1)^(jf + (1-ma)//2) * (1+(-1)^(li+lf))
                            if signQ == 0; continue;end
                            m3j =  wigner3j(Float64,jf//2,lam,ji//2,-ma//2,0,mc//2)
                            coeff = signQ * m3j *Q3j * pfac * ec_p * c_mom_Q * ex_r2
                        end
                        @inbounds for bi = 1:lblock
                            bf = bi # M-conservation
                            TF=[true]; ret=[1,-1]; ridx=[-1,-1,-1]
                            idim = tdims[bi]
                            l_Nn = length(nbits[bi])
                            @inbounds for (Npi,pPhi) in enumerate(pbits[bi])
                                TF_connectable_1(pPhi,bit_a,bit_c,TF)        
                                if TF[1]==false; continue;end
                                calc_phase_1!(pPhi,bit_a,bit_c,ret)
                                bisearch!(pbits[bf],ret[2],ridx)        
                                if ridx[1] == 0;continue;end
                                Npf = ridx[1]
                                if mode =="M1"
                                    push!(pjumps_mu,
                                          op_M1(bi,Npi,Npf,ret[1]==-1,
                                                coeff .*ifelse(Npi==Npf,0.5,1.0)))
                                elseif mode =="E2"
                                    push!(pjumps_Q,
                                          op_E2(bi,Npi,Npf,ret[1]==-1,
                                                coeff .*ifelse(Npi==Npf,0.5,1.0)))
                                end
                            end
                        end
                    end                    
                end
            end
        end
    end
    vec_n_ani = [ false for i=1:length(mstates_n)]
    vec_n_cre = [ false for i=1:length(mstates_n)]
    for i=1:ln
        ni,li,ji,tzi = n_sps[i]
        for f = i:ln
            nf,lf,jf,tzf = n_sps[f]
            ex_r2 = radius_power(2,ni,li,nf,lf,bpar)
            l6j = wigner6j(Float64,lf,li,1,ji//2,jf//2,1//2)
            s6j = wigner6j(Float64,1//2,1//2,1,ji//2,jf//2,lf)
            Q3j = wigner3j(Float64,ji//2,2,jf//2,1//2,0,-1//2)
            jc,mc_s,mc_idxs = possible_mz(n_sps[i],mstates_n)
            ja,ma_s,ma_idxs = possible_mz(n_sps[f],mstates_n)
            for mode in ["M1", "E2"]
                if mode == "M1"
                    if (ni !=nf || li != lf)
                        continue
                    end
                end               
                @inbounds for (ic,mc) in enumerate(mc_s)
                    vec_n_ani .= false
                    vec_n_ani[mc_idxs[ic]] = true
                    bit_c = bitarr_to_int(vec_n_ani)
                    if bit_c == 0; continue; end                
                    @inbounds for (ia,ma) in enumerate(ma_s)
                        if ma != mc; continue;end #M conservation
                        #if a == c && ma > mc; continue;end
                        vec_n_cre .= false
                        vec_n_cre[ma_idxs[ia]] = true
                        bit_a = bitarr_to_int(vec_n_cre)
                        if bit_a == 0; continue; end
                        
                        pfac = sqrt((ji+1)*(jf+1))
                        sign1 = (-1)^((jf-ma)/2+lf+(ji+3)/2)
                        sign2 = (-1)^((jf-ma)/2+lf+(jf+3)/2)
                        if mode =="M1" 
                            fac1 = sign1 * pfac .* sqrt(lf*(lf+1)*(2*lf+1)) * l6j 
                            fac2 = sign2 * pfac * sqrt(1.5) * s6j
                            lam = 1 
                            m3j =  wigner3j(Float64,jf//2,lam,ji//2,-ma//2,0,mc//2)
                            coeff = m3j * (fac1*gfactors[2]+fac2*gfactors[4]) * c_mom_mu
                        elseif mode=="E2"
                            lam = 2
                            signQ = (-1)^(jf + (1-ma)//2) * (1+(-1)^(li+lf))
                            if signQ == 0; continue;end
                            m3j =  wigner3j(Float64,jf//2,lam,ji//2,-ma//2,0,mc//2)
                            coeff = signQ * m3j *Q3j * pfac * ec_n * c_mom_Q * ex_r2
                            #println("fE2 $coeff sign $signQ 3j $m3j jf $jf ji $ji ma $ma mc $mc")
                            if coeff == 0.0; continue;end
                        end
                        @inbounds for bi = 1:lblock
                            bf = bi # M-conservation
                            TF=[true]; ret=[1,-1]; ridx=[-1,-1,-1]
                            nbit = nbits[bi]
                            idim = tdims[bi]; l_Nn = length(nbit)
                            @inbounds for (Nni,nPhi) in enumerate(nbit)
                                TF_connectable_1(nPhi,bit_a,bit_c,TF)
                                if TF[1]==false; continue;end
                                calc_phase_1!(nPhi,bit_a,bit_c,ret)
                                bisearch!(nbits[bf],ret[2],ridx)
                                if ridx[1] == 0;continue;end
                                Nnf = ridx[1]
                                if mode =="M1"
                                    push!(njumps_mu,
                                          op_M1(bi,Nni,Nnf,ret[1]==-1,
                                                coeff .*ifelse(Nni==Nnf,0.5,1.0)))
                                elseif mode =="E2"
                                    push!(njumps_Q,
                                          op_E2(bi,Nni,Nnf,ret[1]==-1,
                                                coeff .*ifelse(Nni==Nnf,0.5,1.0)))
                                end
                            end
                        end
                    end                    
                end
            end
        end
    end
    return pjumps_mu,njumps_mu, pjumps_Q,njumps_Q             
end

function op_obtr_diag(p_oprts,n_oprts,pbits,nbits,tdims,wf,twf)
    @inbounds for tmp in p_oprts
        bi = tmp.bi; Npi = tmp.i; Npf = tmp.f
        coeff = ifelse(tmp.phase,-tmp.fac,tmp.fac)
        pbit = pbits[bi]; nbit = nbits[bi]
        idim = tdims[bi]
        l_Nn = length(nbit)
        l_Np = length(pbit)
        tMi = idim - l_Nn + Npi *l_Nn
        tMf = idim - l_Nn + Npf *l_Nn
        @inbounds for nidx = 1:l_Nn
            Mi = tMi + nidx
            Mf = tMf + nidx
            twf[Mf] += coeff .* wf[Mi]
            twf[Mi] += coeff .* wf[Mf]
        end
    end
    for tmp in n_oprts
        bi = tmp.bi
        Nni = tmp.i; Nnf = tmp.f
        phase = tmp.phase
        coeff = ifelse(phase,-tmp.fac,tmp.fac)
        pbit = pbits[bi]; nbit = nbits[bi]
        idim = tdims[bi]
        l_Nn = length(nbit)
        l_Np = length(pbit)
        tMi = idim - l_Nn + Nni
        tMf = idim - l_Nn + Nnf
        @inbounds for pidx = 1:l_Np
            Mi = tMi + pidx * l_Nn 
            Mf = tMf + pidx * l_Nn 
            twf[Mf] += coeff .* wf[Mi]
            twf[Mi] += coeff .* wf[Mf]
        end
    end
    return nothing
end


function prep_obtr_nondiag(jl2,jr2,p_sps,n_sps,
                           mstates_p,mstates_n,
                           jocc_p,jocc_n,                           
                           tdims_l,pbits_l,nbits_l,
                           tdims_r,pbits_r,nbits_r,
                           p_bifs,n_bifs,omu,bpar,
                           gfactors,eff_charge)
    lp=length(p_sps); ln=length(n_sps);
    lblock_l = length(nbits_l)
    lblock_r = length(nbits_r)
    ec_p,ec_n = eff_charge
    c_mom_mu = sqrt(3.0/(4.0*pi))
    c_mom_Q = sqrt(5.0/(16.0*pi))

    ## M1 (mu; magnetic dipole moment)
    pjumps_mu = op_M1_nd[ ] 
    njumps_mu = op_M1_nd[ ]
    ## E2 (Q; electric quadropole moment)
    pjumps_Q = op_E2_nd[ ]
    njumps_Q = op_E2_nd[ ]

    vec_p_ani = [ false for i=1:length(mstates_p)]
    vec_p_cre = [ false for i=1:length(mstates_p)]
    ex_r2 = 0.0 ; coeff= 0.0
    for i=1:lp
        ni,li,ji,tzi = p_sps[i]
        for f = 1:lp
            nf,lf,jf,tzf = p_sps[f]
            ex_r2 = radius_power(2,ni,li,nf,lf,bpar)
            l6j = wigner6j(Float64,lf,li,1,ji//2,jf//2,1//2)
            s6j = wigner6j(Float64,1//2,1//2,1,ji//2,jf//2,lf)
            Q3j = wigner3j(Float64,ji//2,2,jf//2,1//2,0,-1//2)
            c,mc_s,mc_idxs = possible_mz(p_sps[i],mstates_p)
            a,ma_s,ma_idxs = possible_mz(p_sps[f],mstates_p)
            for mode in ["M1", "E2"]
                if mode == "M1" # slection rule
                    if (ni !=nf || li != lf); continue;end
                    if div(jl2-jr2,2) > 1;continue;end
                elseif mode == "E2"
                    if div(jl2-jr2,2) > 2; continue;end
                end
                @inbounds for (ic,mc) in enumerate(mc_s)
                    vec_p_ani .= false
                    vec_p_ani[mc_idxs[ic]] = true
                    bit_c = bitarr_to_int(vec_p_ani)
                    @inbounds for (ia,ma) in enumerate(ma_s)
                        if div(ma - mc,2) != omu ; continue;end #M conservation
                        vec_p_cre .= false
                        vec_p_cre[ma_idxs[ia]] = true
                        bit_a = bitarr_to_int(vec_p_cre) 
                        pfac = sqrt((ji+1)*(jf+1))
                        if mode =="M1"
                            lam = 1
                            sign1 = (-1)^((jf-ma)/2+lf+(ji+3)/2)
                            sign2 = (-1)^((jf-ma)/2+lf+(jf+3)/2)
                            fac1 = sign1 * l6j .* sqrt(lf*(lf+1)*(2*lf+1))
                            fac2 = sign2 * s6j .* sqrt(1.5)
                            m3j = pfac * wigner3j(jf//2,lam,ji//2,-ma//2,omu,mc//2) 
                            coeff = m3j * (fac1*gfactors[1]+fac2*gfactors[3]) *c_mom_mu
                        elseif mode =="E2"
                            lam = 2
                            signQ = (-1)^(jf + (1-ma)//2) * (1+(-1)^(li+lf))
                            if signQ == 0; continue;end
                            m3j =  wigner3j(Float64,jf//2,lam,ji//2,-ma//2,omu,mc//2)
                            coeff = signQ * m3j *Q3j * pfac * ec_p * c_mom_Q * ex_r2
                        end
                        @inbounds for (bl,br) in p_bifs
                            TF=[true]; ret=[1,-1]; ridx=[-1,-1,-1]
                            @inbounds for (Npi,pPhi) in enumerate(pbits_r[br])
                                TF_connectable_1(pPhi,bit_a,bit_c,TF)        
                                if TF[1]==false; continue;end
                                calc_phase_1!(pPhi,bit_a,bit_c,ret) 
                                bisearch!(pbits_l[bl],ret[2],ridx)        
                                if ridx[1] == 0;continue;end
                                Npf = ridx[1]
                                if mode =="M1"
                                    push!(pjumps_mu,
                                          op_M1_nd(bl,Npf,br,Npi,ret[1]==-1,coeff))
                                elseif mode =="E2"
                                    push!(pjumps_Q,
                                          op_E2_nd(bl,Npf,br,Npi,ret[1]==-1,coeff))
                                end
                            end
                        end
                    end                    
                end
            end
        end
    end
    vec_n_ani = [ false for i=1:length(mstates_n)]
    vec_n_cre = [ false for i=1:length(mstates_n)]
    for i=1:ln
        ni,li,ji,tzi = n_sps[i]
        for f = 1:ln
            nf,lf,jf,tzf = n_sps[f]
            ex_r2 = radius_power(2,ni,li,nf,lf,bpar)
            l6j = wigner6j(Float64,lf,li,1,ji//2,jf//2,1//2)
            s6j = wigner6j(Float64,1//2,1//2,1,ji//2,jf//2,lf)
            Q3j = wigner3j(Float64,ji//2,2,jf//2,1//2,0,-1//2)
            jc,mc_s,mc_idxs = possible_mz(n_sps[i],mstates_n)
            ja,ma_s,ma_idxs = possible_mz(n_sps[f],mstates_n)
            for mode in ["M1", "E2"]
                if mode == "M1"
                    if (ni !=nf || li != lf);continue;end
                    if div(jl2-jr2,2) > 1;continue;end
                elseif mode == "E2"
                    if div(jl2-jr2,2) > 2; continue;end
                end               
                @inbounds for (ic,mc) in enumerate(mc_s)
                    vec_n_ani .= false
                    vec_n_ani[mc_idxs[ic]] = true
                    bit_c = bitarr_to_int(vec_n_ani)
                    if bit_c == 0; continue; end                
                    @inbounds for (ia,ma) in enumerate(ma_s)
                        if div(ma -mc,2) != omu; continue;end #M conservation
                        vec_n_cre .= false
                        vec_n_cre[ma_idxs[ia]] = true
                        bit_a = bitarr_to_int(vec_n_cre)
                        if bit_a == 0; continue;end                        
                        pfac = sqrt((ji+1)*(jf+1))
                        sign1 = (-1)^((jf-ma)/2+lf+(ji+3)/2)
                        sign2 = (-1)^((jf-ma)/2+lf+(jf+3)/2)
                        if mode =="M1" 
                            fac1 = sign1 * pfac .* sqrt(lf*(lf+1)*(2*lf+1)) * l6j 
                            fac2 = sign2 * pfac * sqrt(1.5) * s6j
                            lam = 1 
                            m3j =  wigner3j(Float64,jf//2,lam,ji//2,-ma//2,omu,mc//2)
                            coeff = m3j * (fac1*gfactors[2]+fac2*gfactors[4]) * c_mom_mu
                        elseif mode=="E2"
                            lam = 2
                            signQ = (-1)^(jf + (1-ma)//2) * (1+(-1)^(li+lf))
                            if signQ == 0; continue;end
                            m3j =  wigner3j(Float64,jf//2,lam,ji//2,-ma//2,omu,mc//2)
                            coeff = signQ * m3j *Q3j * pfac * ec_n * c_mom_Q * ex_r2
                        end
                        @inbounds for (bl,br) in n_bifs
                            TF=[true]; ret=[1,-1]; ridx=[-1,-1,-1]
                            @inbounds for (Nni,nPhi) in enumerate(nbits_r[br])
                                TF_connectable_1(nPhi,bit_a,bit_c,TF)
                                if TF[1]==false; continue;end
                                calc_phase_1!(nPhi,bit_a,bit_c,ret)
                                bisearch!(nbits_l[bl],ret[2],ridx)
                                if ridx[1] == 0;continue;end
                                Nnf = ridx[1]
                                if mode =="M1"
                                    push!(njumps_mu,
                                          op_M1_nd(bl,Nnf,br,Nni,ret[1]==-1,coeff))
                                elseif mode =="E2"
                                    push!(njumps_Q,
                                          op_E2_nd(bl,Nnf,br,Nni,ret[1]==-1,coeff))
                                end
                            end
                        end
                    end                    
                end
            end
        end
    end
    return pjumps_mu,njumps_mu,pjumps_Q,njumps_Q             
end

function op_obtr_nondiag(p_oprts,n_oprts,
                         pbits_l,nbits_l,tdims_l,
                         pbits_r,nbits_r,tdims_r,
                         wf,twf)
    @inbounds for tmp in p_oprts
        bi = tmp.br; Npi = tmp.Nr
        bf = tmp.bl; Npf = tmp.Nl
        coeff = ifelse(tmp.phase,-tmp.fac,tmp.fac)
        nbit_r = nbits_r[bi]; nbit_l = nbits_l[bf]
        idim = tdims_r[bi]
        fdim = tdims_l[bf]
        l_Nn_r = length(nbit_r)
        l_Nn_l = length(nbit_l)
        tMi = idim - l_Nn_r + Npi *l_Nn_r
        tMf = fdim - l_Nn_l + Npf *l_Nn_l
        @inbounds for nidx_r = 1:l_Nn_r
            Mi = tMi + nidx_r
            n_r = nbit_r[nidx_r]
            @inbounds for nidx_l = 1:l_Nn_l
                n_l = nbit_l[nidx_l]
                if n_r != n_l; continue;end
                Mf = tMf + nidx_l
                twf[Mf] += coeff .* wf[Mi]
            end
        end
    end
    @inbounds for tmp in n_oprts
        bi = tmp.br; Nni = tmp.Nr
        bf = tmp.bl; Nnf = tmp.Nl
        coeff = ifelse(tmp.phase,-tmp.fac,tmp.fac)
        pbit_r = pbits_r[bi]; pbit_l = pbits_l[bf]
        nbit_r = nbits_r[bi]; nbit_l = nbits_l[bf]
        idim = tdims_r[bi]
        fdim = tdims_l[bf]
        l_Nn_r = length(nbit_r)
        l_Nn_l = length(nbit_l)
        l_Np_r = length(pbit_r)
        l_Np_l = length(pbit_l)
        tMi = idim - l_Nn_r + Nni
        tMf = fdim - l_Nn_l + Nnf
        @inbounds for pidx_r = 1:l_Np_r
            Mi = tMi + pidx_r * l_Nn_r
            p_r = pbit_r[pidx_r]
            @inbounds for pidx_l = 1:l_Np_l
                p_l = pbit_l[pidx_l]
                if p_r != p_l; continue;end
                Mf = tMf + pidx_l * l_Nn_l
                twf[Mf] += coeff .* wf[Mi]
            end
        end
    end
    return nothing
end

function calc_bifs(Mps_l,Mns_l,Mtot_l,Mps_r,Mns_r,Mtot_r)
    # for right wav
    blocks_r = [ [0,0] ]; deleteat!(blocks_r,1)
    for Mp in Mps_r
        for Mn in Mns_r
            if Mp + Mn == Mtot_r
                push!(blocks_r,[Mp,Mn])
            end
        end
    end
    blocks_l = [ [0,0] ]; deleteat!(blocks_l,1)
    for Mp in Mps_l
        for Mn in Mns_l
            if Mp + Mn == Mtot_l
                push!(blocks_l,[Mp,Mn])
            end
        end
    end
    p_bifs = [ [0,0] ]; deleteat!(p_bifs,1)
    n_bifs = [ [0,0] ]; deleteat!(n_bifs,1)

    for (ir,bl_r) in enumerate(blocks_r)
        Mpr,Mnr = bl_r
        for (il,bl_l) in enumerate(blocks_l)
            Mpl,Mnl = bl_l
            if Mnr == Mnl
                push!(p_bifs,[il,ir])
            end
            if Mpr == Mpl
                push!(n_bifs,[il,ir])
            end
        end
    end
    return p_bifs,n_bifs
end

"""
    transit_main(sntf,target_nuc,jl2,jr2,in_wfs;
                 num_ev_l=100,num_ev_r=100,q=1,is_block=false,is_show=true,
                 calc_EM=true,gfactors=[1.0,0.0,5.586,-3.826],eff_charge=[1.5,0.5])

Calculate the M1&E2 transitions for two wavefunctions

# Arguments
- `sntf`:path to the interaction file
- `target_nuc`:target nucleus in string e.g., "Si28"
- `jl2`:J*2 for the left w.f.
- `jr2`:J*2 for the right w.f.
- `in_wfs`:["path to left wf","path to right wf"]

# Optional arguments
- `num_ev_l=100`:upper limit of the number of eigenvectors for the left w.f.
- `num_ev_r=100`:upper limit of the number of eigenvectors for the right w.f.
- `is_show=true`:to display the TimerOutput
- `gfactors=[1.0,0.0,5.586,-3.826]`:g factors [glp,gln,gsp,gsn]
- `eff_charge=[1.5,0.5]`:effective charges [ep,en]

# Optional arguments (not used, but for future developments)
- `q=1`:block size
- `is_block=false`:to use Block algorithm
"""
function transit_main(sntf,target_nuc,jl2,jr2,in_wfs;
                      num_ev_l=100,num_ev_r=100,
                      q=1,is_block=false,is_show=true,
                      calc_EM=true,
                      gfactors=[1.0,0.0,5.586,-3.826],
                      eff_charge=[1.5,0.5])
    if length(in_wfs) == 0;println("input wfs must be specified");exit();end
    
    to = TimerOutput()
    Anum = parse(Int64, match(reg,target_nuc).match)
    lp,ln,cp,cn,massop,Aref,pow,p_sps,n_sps,SPEs,olabels,oTBMEs,labels,TBMEs = readsnt(sntf,Anum)
    hw, bpar = init_ho_by_mass(Anum,1) # mass formula 
    if 16 <= Anum <= 40
        hw, bpar = init_ho_by_mass(Anum,2) # 2: mass formula for sd-shell
    end

    target_el = replace.(target_nuc, string(Anum)=>"")
    Z,N,vp,vn = getZNA(target_el,Anum,cp,cn)
    mstates_p, mstates_n,mz_p,mz_n = def_mstates(p_sps,n_sps)

    widx = [1,2]
    if jl2 < jr2
        tr = jl2; tl = jr2; jl2 = tl; jr2 = tr; widx=[2,1]
    end
    
    Mtot_l = jl2; Mtot_r = jr2
    pbits_l,nbits_l,jocc_p,jocc_n,Mps_l,Mns_l,tdims_l = occ(p_sps,mstates_p,mz_p,vp,
                                                            n_sps,mstates_n,mz_n,vn,
                                                            Mtot_l)
    pbits_r,nbits_r,jocc_p,jocc_n,Mps_r,Mns_r,tdims_r = occ(p_sps,mstates_p,mz_p,vp,
                                                            n_sps,mstates_n,mz_n,vn,
                                                            Mtot_r)
    mdim_l = tdims_l[end];   mdim_r = tdims_r[end]
    println("wfl: ",in_wfs[widx[1]])
    if in_wfs[widx[1]] != in_wfs[widx[2]]; println("wfr: ",in_wfs[widx[2]]);end
    println("mdims J2=$jl2 $mdim_l J2=$jr2 $mdim_r Mtots=($Mtot_l, $Mtot_r)")
    wfs_l,jjs_l = read_wav(in_wfs[widx[1]],mdim_l,1;all=true)
    wfs_r,jjs_r = read_wav(in_wfs[widx[2]],mdim_r,1;all=true)
    lblock_l=length(pbits_l); lblock_r=length(pbits_r)
    omu = div(jl2-jr2,2)
    if omu > 2; prinln("Only the M1&E2 are supported now");return nothing;end
   
    p_bifs,n_bifs = calc_bifs(Mps_l,Mns_l,Mtot_l,Mps_r,Mns_r,Mtot_r)
    # operation of EM transition operators O|Phi>
    @timeit to "prep" begin
        op_p_mu,op_n_mu,op_p_Q,op_n_Q = prep_obtr_nondiag(jl2,jr2,p_sps,n_sps,
                                                          mstates_p,mstates_n,
                                                          jocc_p,jocc_n,
                                                          tdims_l,pbits_l,nbits_l,
                                                          tdims_r,pbits_r,nbits_r,
                                                          p_bifs,n_bifs,omu,bpar,
                                                          gfactors,eff_charge)
    end

    ## eval expectation values
    c_mom_mu = sqrt(4.0*pi/3.0)
    c_mom_Q  = sqrt(16.0*pi/5.0)         
    
    if tri_relation(jl2//2,jr2//2,1)
        @timeit to "M1" begin
        #println("M1 available")
        twf = zeros(Float64,mdim_l)
        lam = 1
        Jf = jl2//2; Ji=jr2//2
        Mf=Mtot_l//2;Mi=Mtot_r//2
        sel = wigner3j(Jf,lam,Ji,-Mf,omu,Mi)
        signmu = (-1)^(Jf-Mf)
        cg = wigner3j(Ji,lam,Jf,-Ji,omu,Jf) 
        for (nth_r,wf_r) in enumerate(wfs_r)
            if sel == 0; continue;end
            if nth_r > num_ev_r; continue;end
            twf .= 0.0
            op_obtr_nondiag(op_p_mu,op_n_mu,                     
                            pbits_l,nbits_l,tdims_l,
                            pbits_r,nbits_r,tdims_r,
                            wf_r,twf)
            if calc_EM 
                for (nth_l,wf_l) in enumerate(wfs_l)
                    if nth_l > num_ev_l; continue;end
                    Mred =  dot(wf_l,twf) /sel
                    tx = ""
                    tx *= "Mred.(M1) ($nth_l,$nth_r) " * @sprintf("%9.3f ",Mred) *"  "
                    B_l = Mred^2 / (2*Ji+1); B_r = Mred^2 / (2*Jf+1)
                    tx *= "B(M1) ($nth_l<=$nth_r) " * @sprintf("%9.3f ",B_l) *" B(M1) ($nth_l=>$nth_r) " * @sprintf("%9.3f ",B_r)
                    if in_wfs[1] == in_wfs[2] && nth_l==nth_r
                        gm = Mred * cg *c_mom_mu * signmu
                        tx *= " <gm>= "* @sprintf("%9.3f ",gm)
                    end
                    println(tx)
                end
            else
                nth_l = nth_r
                wf_l= wfs_l[nth_l]
                Mred =  dot(wf_l,twf) /sel
                gm = Mred * cg *c_mom_mu * signmu
                tx = " <gm>= "* @sprintf("%9.3f ",gm) 
                println(tx)
            end
        end
        end
    end
    if tri_relation(jl2//2,jr2//2,2)
        #println("E2 available")
        @timeit to "E2" begin
        twf = zeros(Float64,mdim_l)
        lam = 2
        Jf = jl2//2; Mf = Mtot_l//2
        Ji = jr2//2; Mi = Mtot_r//2        
        sel = wigner3j(Jf,lam,Ji,-Mf,omu,Mi)
        signQ = (-1)^(Jf-Mf)
        cg = wigner3j(Ji,lam,Jf,-Ji,omu,Jf) 
        for (nth_r,wf_r) in enumerate(wfs_r)
            if sel==0;continue;end
            if nth_r > num_ev_r; continue;end
            twf .= 0.0
            op_obtr_nondiag(op_p_Q,op_n_Q,   
                            pbits_l,nbits_l,tdims_l,
                            pbits_r,nbits_r,tdims_r,
                            wf_r,twf)
            if calc_EM 
                for (nth_l,wf_l) in enumerate(wfs_l)
                    if nth_l > num_ev_l; continue;end
                    Mred =  dot(wf_l,twf) /sel
                    tx = ""
                    tx *= "Mred.(E2) ($nth_l,$nth_r) " * @sprintf("%9.3f ",Mred) *"  "
                    B_l = Mred^2 / (2*Ji+1); B_r = Mred^2 / (2*Jf+1)
                    tx *= "B(E2) ($nth_l<=$nth_r) " * @sprintf("%9.3f ",B_l) *" B(E2) ($nth_l=>$nth_r) " * @sprintf("%9.3f ",B_r)
                    if in_wfs[1] == in_wfs[2] && nth_l==nth_r
                        Q = Mred * cg *c_mom_Q * signQ
                        tx *= " <eQ>= "* @sprintf("%9.3f ",Q)
                    end
                    println(tx)
                end
            else
                wf_l = wfs_l[nth_r]
                Mred =  dot(wf_l,twf) /sel
                Q = Mred * cg *c_mom_Q * signQ
                tx = " <eQ>= "* @sprintf("%9.3f ",Q)
                println(tx)
            end
        end
        end
    end
    if is_show
        show(to, allocations = true,compact = false)
    end
    println("")
    return nothing
end

function eval_moment(Mtot,Rvecs,totJs,
                     p_sps,n_sps,
                     mstates_p,mstates_n,tdims,
                     jocc_p,jocc_n,pbits,nbits,bpar,
                     gfactors,eff_charge)
    mdim = tdims[end]
    c_mom_mu = sqrt(4.0*pi/3.0)
    c_mom_Q = sqrt(16.0*pi/5.0)        
    op_p_mu,op_n_mu,op_p_Q,op_n_Q = prep_obtr_diag(p_sps,n_sps,
                                                   mstates_p,mstates_n,
                                                   tdims,jocc_p,jocc_n,
                                                   pbits,nbits,bpar,
                                                   gfactors,eff_charge)
    tx_mom = ""
    twf1 = zeros(Float64,mdim)
    twf2 = zeros(Float64,mdim)
    Mi = Mtot//2; Mf = Mtot//2
    for (nth,wf) in enumerate(Rvecs)
        twf1 .= 0.0 ; twf2 .= 0.0
        # for M1
        p_oprts = op_p_mu;  n_oprts = op_n_mu
        op_obtr_diag(p_oprts,n_oprts,pbits,nbits,tdims,wf,twf1)
        # for E2
        p_oprts = op_p_Q;  n_oprts = op_n_Q
        op_obtr_diag(p_oprts,n_oprts,pbits,nbits,tdims,wf,twf2)
        
        Ji = Int(2*totJs[nth]) // 2; 
        for np = nth:length(Rvecs)
            tx = ""
            wf_l = Rvecs[np]
            Jf = Int(2*totJs[np])//2
            lam = 1
            omu = div(Mf-Mi,2)
            sel = wigner3j(Float64,Jf,lam,Ji,-Mf,omu,Mi)
            if sel != 0
                signmu = (-1)^(Jf-Mf)
                cg = wigner3j(Float64,Ji,lam,Jf,-Ji,omu,Jf) / sel
                if np ==nth
                    gm = signmu * dot(twf1,wf) *cg  * c_mom_mu
                    tx *= "mu (n=$nth) " * @sprintf("%12.5f ",gm) *"\t"
                else
                    Mred =  dot(wf_l,twf1) /sel
                    B_l = Mred^2 / (2*Ji+1); B_r = Mred^2 / (2*Jf+1)
                    #tx *= "Mred.(M1) ($np,$nth) " * @sprintf("%9.3f ",Mred) *"  "      
                    #tx *= "B(M1) ($np<=$nth) " * @sprintf("%9.3f ",B_l)
                    #tx *= " B(M1) ($np=>$nth) " * @sprintf("%9.3f ",B_r)
                end
            end
            lam = 2
            sel =  wigner3j(Float64,Jf,lam,Ji,-Mf,omu,Mi)
            if sel != 0
                signQ = ( (-1)^(Jf-Mf))
                cg = wigner3j(Float64,Ji,lam,Jf,-Ji,omu,Jf) / sel
                if np == nth                        
                    Q = dot(twf2,wf) * cg * c_mom_Q * signQ
                    tx *= "Q (n=$nth) " * @sprintf("%12.5f ",Q)
                else
                    Mred = dot(wf_l,twf2) /sel
                    B_l = Mred^2 / (2*Ji+1); B_r = Mred^2 / (2*Jf+1)
                    #tx *= "Mred.(E2) ($np,$nth) " * @sprintf("%9.3f ",Mred) *"  "
                    #tx *= "B(E2) ($np<=$nth) " * @sprintf("%9.3f ",B_l)
                    #tx *= " B(E2) ($np=>$nth) " * @sprintf("%9.3f ",B_r) 
                end
            end
            if tx !="";tx_mom *= tx*"\n";end
        end
    end
    return tx_mom
end
