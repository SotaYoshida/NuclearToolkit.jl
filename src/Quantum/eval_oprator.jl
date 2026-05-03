function add_HCM_to_Hmat!(Hmat, HCMmat)
    for (nkey, val) in HCMmat
        if haskey(Hmat, nkey)
            Hmat[nkey] += val
        else
            Hmat[nkey] = val
        end
    end
    return nothing
end

function eval_HCMexpec(beta_cm, evecs, mdim, HCMmat, Ecms, n_eigen, to)
    v = zeros(Float64, mdim)
    Hv = zeros(Float64, mdim)
    partials = zeros(Float64, mdim, Threads.maxthreadid())
    keysvec = collect(keys(HCMmat))
    for i in 1:n_eigen
        v .= evecs[:, i]
        operate_H_on_vec!(Hv, HCMmat, v, partials, keysvec, to)
        tEcm = dot(v, Hv)
        if isnan(tEcm)
            println("Ecm computation resulted in NaN for state $i")
            println("norm(v): $(norm(v)) norm(Hv): $(norm(Hv))")
            println("v: $v")
        end
        Ecms[i] = tEcm
    end
    print_vec("Ecm (MeV):     ", Ecms ./ beta_cm)
    return nothing
end

function construct_HCMmat(beta_cm, mdim, all_bitint_prod, p_msps, n_msps, vZ, vN, 
                          dict_m2j, int_shift, dict_CGs, dict_HCM, verbose, to;
                          Hrank=2)
    HCMmat = Dict{UInt, Float64}( )
    partial_dict = [ Dict{UInt, Float64}( ) for _ in 1:Threads.maxthreadid() ]
    possible_ijkls_mthread = [ Vector{Vector{Int}}() for _ in 1:Threads.maxthreadid() ]
    @timeit to "thread loop" Threads.@threads for iter in 1:mdim^2
        possible_ijkls_worker = possible_ijkls_mthread[Threads.threadid()]
        idx_bra = div(iter-1, mdim) + 1
        idx_ket = iter - (idx_bra-1) * mdim
        if idx_bra > idx_ket
            continue # skip lower triangle
        end

        p_bitint_bra, n_bitint_bra = all_bitint_prod[idx_bra]
        p_bitint_ket, n_bitint_ket = all_bitint_prod[idx_ket]

        p_ham_dist = count_ones( p_bitint_bra ⊻ p_bitint_ket )
        n_ham_dist = count_ones( n_bitint_bra ⊻ n_bitint_ket )

        sum_ham_dist = p_ham_dist + n_ham_dist
        if sum_ham_dist > Hrank * 2 # bra and ket cannot be connected by the Hamiltonian
            continue 
        end

        tid = Threads.threadid()
        # # Calculate the Hamiltonian matrix element
        Hij = 0.0

        #@timeit to "diag-1b" 
        if sum_ham_dist == 0
            Hij += beta_cm * eval_HCM_1b(p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket,
                                         p_msps, n_msps, dict_m2j, int_shift, dict_CGs)                
        end

        #@timeit to "diag-2b" 
        if sum_ham_dist == 0 
            Hij += beta_cm * eval_HCM_2b_diag(p_bitint_bra, n_bitint_bra, p_msps, n_msps, dict_CGs, dict_HCM, int_shift, verbose)            
        end

        #@timeit to "nondiag-2b" 
        if sum_ham_dist == 2 || sum_ham_dist == 4 
            Hij += beta_cm * eval_HCM_2b_nond(p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket,
                                              p_ham_dist, n_ham_dist,
                                              p_msps, n_msps, possible_ijkls_worker, dict_CGs, dict_HCM, int_shift, to, verbose)
        end

        #@timeit to "assign" 
        partial_dict[tid][hash_2ints(idx_bra, idx_ket)] = Hij
    end    

    # Combine the partial results from all threads
    @timeit to "collect" for tid in 1:Threads.maxthreadid()
        for (nkey, val) in partial_dict[tid]
            HCMmat[nkey] = val
        end
    end

    if verbose >=1
        println("HCM matrix")
        for (nkey, val) in HCMmat
            idx_bra, idx_ket = unhash_2ints(nkey)
            bra_bitint, ket_bitint = all_bitint_prod[idx_bra], all_bitint_prod[idx_ket]
            p_bra, n_bra = bra_bitint
            p_ket, n_ket = ket_bitint
            bra_bitstr = int2bitstr(p_bra, length(p_msps)) * " ⊗ " * int2bitstr(n_bra, length(n_msps))
            ket_bitstr = int2bitstr(p_ket, length(p_msps)) * " ⊗ " * int2bitstr(n_ket, length(n_msps))
            println("Hcm[$(idx_bra),$(idx_ket)] = $(HCMmat[nkey]) for bra $(bra_bitstr) ket $(ket_bitstr)")
        end
    end
    return HCMmat
end

"""
<sps(n,l,j,tz)|| Hcm_1b || sps(n,l,j,tz)> diagonal part of the HCM 1-body term
"""
function hcm_1body(nljt1, nljt2, idxs_core)
    # Tuple-like case 
    n1, l1, j1, t1 = nljt1
    n2, l2, j2, t2 = nljt2
    r = 0.0
    if n1 != n2 || l1 != l2 || j1 != j2 || t1 !=t2
        return r
    end
    r = (2*n1 + l1) 
    # For NCSM, ended here.
    
    #* sqrt(j1+1) # 
    for idx in idxs_core
        state = msps[idx]
        nc = mstate.n; lc = mstate.l; jc = mstate.l; tc = mstate.tz
        Jmin = div(abs(jc-j1),2); Jmax = div(abs(jc+j1),2)
        for Jtot in Jmin:Jmax
            r += hcm_2body(nc, lc, jc, tc, n1, l1, j1, t1,
                           nc, lc, jc, tc, n2, l2, j2, t2, Jtot) / sqrt(j1+1)
        end
    end

   return r 
end

"""
<sps_1, sps_2|| Hcm_2b ||sps_3, sps_4>
"""
function hcm_2body(n1, l1, j1, t1, n2, l2, j2, t2,
                   n3, l3, j3, t3, n4, l4, j4, t4, JJ; verbose=false)
    r = 0.0
    if (t1+t2 != t3+t4) || ((l1+l2)%2 != (l3+l4)%2 ) #isospin/parity check
        return r
    end
    r = (-1)^(div(j2+j3,2)+JJ) * wigner6j(Float64, j1/2, j2/2, JJ, j4/2, j3/2, 1) *(
          - nabla_j(n1,l1,j1,n3,l3,j3)*nabla_j(n2,l2,j2,n4,l4,j4) 
          + radius_j(n1,l1,j1,n3,l3,j3)*radius_j(n2,l2,j2,n4,l4,j4))
    if t3 == t4 
        r -= (-1)^(div(j3+j4,2)-JJ) * (-1)^(div(j2+j4,2)+JJ) * wigner6j(Float64, j1/2, j2/2, JJ, j3/2, j4/2, 1) * (
            - nabla_j(n1,l1,j1,n4,l4,j4)*nabla_j(n2,l2,j2,n3,l3,j3) 
            + radius_j(n1,l1,j1,n4,l4,j4)*radius_j(n2,l2,j2,n3,l3,j3))
    end
    return r
end

"""
<nlj|| ∇*b || nlj>

Note that j1 and j2 are doubled.
"""
function nabla_j(n1, l1, j1, n2, l2, j2)
    r = red_nabla_l(n1, l1, n2, l2)
    r *= (-1)^( l1+div(j2+1, 2)) * sqrt( (j1+1)*(j2+1)) * wigner6j(Float64, j1/2, 1, j2/2, l2, 1/2, l1)
    return r
end

"""
<nlj|| r || nlj>
"""
function radius_j(n1, l1, j1, n2, l2, j2)
    r = 0.0
    if (n1==n2 && l1==l2-1) 
       r = -sqrt(l2*(n2 + l2 + 0.5))
    elseif (n1==n2-1 && l1==l2+1) 
       r = -sqrt((l2+1.0)*n2)
    elseif (n1==n2+1 && l1 == l2-1) 
       r =  sqrt(l2*(n2 + 1.0))
    elseif (n1==n2 && l1==l2+1)
       r =  sqrt((l2 + 1.0)*(n2 + l2 + 1.5))
    else 
       r = 0.0
    end
    r *= (-1)^( l1+div(j2+1, 2)) * sqrt( (j1+1)*(j2+1) ) * wigner6j(Float64, j1/2, 1, j2/2, l2, 1/2, l1)
    return r
end

function eval_HCM_1b(p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket,
                     p_msps, n_msps, dict_m2j, int_shift, dict_CGs)::Float64
    # One-body term is diagonal
    H_1b = 0.0
    for idx_m in 1:length(p_msps)
        idx_j = dict_m2j[idx_m]
        onehot = Int128(1) << (idx_m-1)
        if (p_bitint_bra & onehot) == onehot == (p_bitint_ket & onehot) 
            sps = p_msps[idx_m]
            n = sps.n; l = sps.l; j = sps.j; tz = sps.tz
            H_1b += hcm_1body((n, l, j, tz), (n, l, j, tz), Int[ ])
        end
    end
    offset = length(p_msps)
    for idx_m in 1:length(n_msps)
        idx_j = dict_m2j[idx_m+offset]
        onehot = Int128(1) << (idx_m-1)
        if (n_bitint_bra & onehot) == onehot == (n_bitint_ket & onehot) 
            sps = n_msps[idx_m]
            n = sps.n; l = sps.l; j = sps.j; tz = sps.tz
            H_1b += hcm_1body((n, l, j, tz), (n, l, j, tz), Int[ ])
        end
    end
    return H_1b 
end


"""
Returns possible m-scheme indices (i, j, k, l) for the given channel and difference bitstrings.

For hamming distance 2, it is a bit tricky since we need to consider spectator indices.
"""
function get_possible_ijkls!(possible_ijkls,
                             channel, p_msps, n_msps,
                             p_diff, n_diff, count_p_diff, count_n_diff,
                             p_bitint_bra, n_bitint_bra,
                             p_bitint_ket, n_bitint_ket
                            )
    empty!(possible_ijkls)
    msps_1 = ifelse(channel == "nn", n_msps, p_msps)
    ket_1 = ifelse(channel == "nn", n_bitint_ket, p_bitint_ket)
    bra_1 = ifelse(channel == "nn", n_bitint_bra, p_bitint_bra)
    diff_1 = ifelse(channel == "nn", n_diff, p_diff)
    crea_1 = diff_1 & bra_1
    anih_1 = diff_1 & ket_1

    msps_2 = ifelse(channel == "pp", p_msps, n_msps)
    diff_2 = ifelse(channel == "pp", p_diff, n_diff)
    bra_2 = ifelse(channel == "pp", p_bitint_bra, n_bitint_bra)
    ket_2 = ifelse(channel == "pp", p_bitint_ket, n_bitint_ket)
    crea_2 = diff_2 & bra_2
    anih_2 = diff_2 & ket_2
    hamm_dist = count_p_diff + count_n_diff

    # First, asign naively. Then, we will count possible permutations later.
    idx_bra_1 = idx_bra_2 = idx_ket_1 = idx_ket_2 = 0
    if hamm_dist == 2
        @assert channel != "pn" "pn channel cannot have hamm dist 2"
        for idx_m in 1:length(msps_1)
            onehot = Int128(1) << (idx_m-1)
            if (crea_1 & onehot) != 0
                idx_bra_1 = idx_m
            end
            if (anih_2 & onehot) != 0
                idx_ket_1 = idx_m
            end
        end
        tf = idx_bra_1 != 0 && idx_ket_1 != 0
        @assert tf "Failed to identify nond. creation/annihilation indices for hamm dist 2. $idx_bra_1 $idx_ket_1"
        for idx_spec = 1:length(msps_2)
            onehot = Int128(1) << (idx_spec-1)
            occupied = (bra_2 & onehot) == onehot == (ket_2 & onehot)
            if !occupied
                continue
            end
            # spectator index must be coupled with both idx_bra_1 and idx_ket_1 via Hcm_2body operators
            n_bra, l_bra = get_nl(msps_1, idx_bra_1)
            n_ket, l_ket= get_nl(msps_1, idx_ket_1)
            n_spec, l_spec = get_nl(msps_1, idx_spec)
            if red_nabla_l(n_bra, l_bra, n_spec, l_spec) == 0.0
                continue
            end
            if red_nabla_l(n_ket, l_ket, n_spec, l_spec) == 0.0
                continue
            end
            # We must consider just whether spectator index can be coupled with both bra and ket indices,
            # and then order them. Exchange terms are handled in `hcm_2body`.
            push!(possible_ijkls, [min(idx_bra_1, idx_spec), max(idx_bra_1, idx_spec),
                                          min(idx_ket_1, idx_spec), max(idx_ket_1, idx_spec)])
        end
    elseif hamm_dist == 4
        hit_count_cre = hit_count_ani = 0
        if channel == "pp" || channel == "nn"
            for idx_m in 1:length(msps_1)
                onehot = Int128(1) << (idx_m-1)
                if (crea_1 & onehot) != 0
                    hit_count_cre += 1
                    if hit_count_cre == 1
                        idx_bra_1 = idx_m
                    elseif hit_count_cre == 2
                        idx_bra_2 = idx_m
                    else
                        error("Too many creation operators found for hamm dist 4")
                    end
                end
                if (anih_1 & onehot) != 0
                    hit_count_ani += 1
                    if hit_count_ani == 1
                        idx_ket_1 = idx_m
                    elseif hit_count_ani == 2
                        idx_ket_2 = idx_m
                    else
                        error("Too many annihilation operators found for hamm dist 4")
                    end
                end
            end
        else
            for idx_m in 1:length(msps_1)
                onehot = Int128(1) << (idx_m-1)
                if (crea_1 & onehot) != 0
                    idx_bra_1 = idx_m
                end
                if (anih_1 & onehot) != 0
                    idx_ket_1 = idx_m
                end
            end
            for idx_m in 1:length(msps_2)
                onehot = Int128(1) << (idx_m-1)
                if (crea_2 & onehot) != 0
                    idx_bra_2 = idx_m
                end
                if (anih_2 & onehot) != 0
                    idx_ket_2 = idx_m
                end
            end
        end
        tf = idx_bra_1 != 0 && idx_bra_2 != 0 && idx_ket_1 != 0 && idx_ket_2 != 0
        @assert tf "Failed to identify nond. creation/annihilation indices for hamm dist 4. $idx_bra_1 $idx_bra_2 $idx_ket_1 $idx_ket_2"
        # If channel == "pn", both proton and neutron parts are unique, so no permutation is needed at the first place.
        # For pp and nn channels, we need to consider exchange terms, but that is already handled in the `hcm_2body` function.
        push!(possible_ijkls, [idx_bra_1, idx_bra_2, idx_ket_1, idx_ket_2])
        return nothing
    else
        error("Invalid hamming distance: $hamm_dist")
    end
    return nothing
end

"""
only a^†_i a^†_j a_j a_i terms (pp/nn) can contribute
"""
function eval_HCM_2b_diag(p_bitint, n_bitint, p_msps, n_msps, dict_CGs, dict_HCM, int_shift, verbose)
    ret = 0.0
    for pn_idx = 1:2
        msps = if pn_idx == 1 p_msps else n_msps end
        bitint = if pn_idx == 1 p_bitint else n_bitint end
        for i = 1:length(msps)
            if (bitint & (Int128(1) << (i-1))) == 0 # i should be occupied
                continue
            end
            n_i = msps[i].n; l_i = msps[i].l; j_i = msps[i].j; jz_i = msps[i].jz; tz_i = msps[i].tz
            key_i = get_nkey4_shift(n_i, l_i, j_i, tz_i)
            for j = i+1:length(msps)
                if (bitint & (Int128(1) << (j-1))) == 0 # j should be occupied
                    continue
                end
                n_j = msps[j].n; l_j = msps[j].l; j_j = msps[j].j; jz_j = msps[j].jz; tz_j = msps[j].tz
                key_j = get_nkey4_shift(n_j, l_j, j_j, tz_j)
                Jmin = div(abs(j_i-j_j), 2)
                Jmax = div(abs(j_i+j_j), 2)
        
                for Jtot in Jmin:Jmax
                    if abs(jz_i + jz_j) > Jtot * 2
                        continue
                    end
                    nkey = (key_i, key_j, key_i, key_j, UInt(Jtot))
                    tmp = get(dict_HCM, nkey, 0.0)
                    cg1 = dict_CGs[get_nkey6_shift(j_i, jz_i, j_j, jz_j, Jtot, div(jz_i+jz_j, 2); int_shift=int_shift)]
                    ret += tmp * cg1^2 
                end
            end
        end
    end
    return ret
end


function get_anticomm_phase(channel, bra1_bit, bra2_bit, ket1_bit, ket2_bit,
                            idx_m_bra_1, idx_m_bra_2,
                            idx_m_ket_1, idx_m_ket_2)
    phase_bra = phase_ket = 1.0
    if channel == "pn"
        for idx_m = 1:idx_m_bra_1-1
            onehot = Int128(1) << (idx_m-1)
            if bra1_bit & onehot == onehot
                phase_bra *= -1.0
            end
        end
        for idx_m in 1:idx_m_ket_1-1
            onehot = Int128(1) << (idx_m-1)
            if ket1_bit & onehot == onehot
                phase_ket *= -1.0
            end
        end
        for idx_m in 1:idx_m_bra_2-1  
            onehot = Int128(1) << (idx_m-1)
            if bra2_bit & onehot == onehot
                phase_bra *= -1.0
            end
        end
        for idx_m in 1:idx_m_ket_2-1
            onehot = Int128(1) << (idx_m-1)
            if ket2_bit & onehot == onehot
                phase_ket *= -1.0
            end
        end
    else
        for idx_m in idx_m_bra_1+1 : idx_m_bra_2-1
            onehot = Int128(1) << (idx_m-1)
            if (bra1_bit & onehot) == onehot == (bra2_bit & onehot)
                phase_bra *= -1.0
            end
        end
        for idx_m in idx_m_ket_1+1 : idx_m_ket_2-1
            onehot = Int128(1) << (idx_m-1)
            if (ket1_bit & onehot) == onehot == (ket2_bit & onehot)
                phase_ket *= -1.0
            end
        end
    end
    return phase_bra, phase_ket
end


function eval_HCM_2b_nond(p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket,
                          count_p_diff, count_n_diff,
                          p_msps, n_msps, 
                          possible_ijkls::Vector{Vector{Int}},
                          dict_CGs::Dict{UInt64, Float64},
                          dict_HCM::Dict{Tuple{UInt, UInt, UInt, UInt, UInt}, Float64},
                          int_shift::Int,
                          to::TimerOutput,
                          verbose::Int):: Float64
    ret = 0.0
    p_diff = p_bitint_bra ⊻ p_bitint_ket
    n_diff = n_bitint_bra ⊻ n_bitint_ket

    channel = "pp"
    if count_p_diff == 4 || (count_p_diff == 2 && count_n_diff == 0)
        channel = "pp"
    elseif count_p_diff == 2 && count_n_diff == 2
        channel = "pn"        
    elseif count_n_diff == 4 || (count_n_diff == 2 && count_p_diff == 0)
        channel = "nn"
    else
        error("Invalid configuration: p $(count_p_diff) n $(count_n_diff)")
    end

    get_possible_ijkls!(possible_ijkls, 
                        channel, p_msps, n_msps, p_diff, n_diff, 
                        count_p_diff, count_n_diff,
                        p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket)

    msps_1 = ifelse(channel == "nn", n_msps, p_msps)
    msps_2 = ifelse(channel == "pp", p_msps, n_msps)
    bra1_bit = ifelse(channel == "nn", n_bitint_bra, p_bitint_bra)
    bra2_bit = ifelse(channel == "pp", p_bitint_bra, n_bitint_bra)
    ket1_bit = ifelse(channel == "nn", n_bitint_ket, p_bitint_ket)
    ket2_bit = ifelse(channel == "pp", p_bitint_ket, n_bitint_ket)

    #@timeit to "loop over possible_ijkls" 
    for ijkl in possible_ijkls
        idx_m_bra_1, idx_m_bra_2, idx_m_ket_1, idx_m_ket_2 = ijkl
        sps_bra_1 = msps_1[idx_m_bra_1]
        sps_bra_2 = msps_2[idx_m_bra_2]
        sps_ket_1 = msps_1[idx_m_ket_1]
        sps_ket_2 = msps_2[idx_m_ket_2]

        n1, l1, j1, j1z, t1 = sps_bra_1.n, sps_bra_1.l, sps_bra_1.j, sps_bra_1.jz, sps_bra_1.tz
        n3, l3, j3, j3z, t3 = sps_bra_2.n, sps_bra_2.l, sps_bra_2.j, sps_bra_2.jz, sps_bra_2.tz
        n2, l2, j2, j2z, t2 = sps_ket_1.n, sps_ket_1.l, sps_ket_1.j, sps_ket_1.jz, sps_ket_1.tz
        n4, l4, j4, j4z, t4 = sps_ket_2.n, sps_ket_2.l, sps_ket_2.j, sps_ket_2.jz, sps_ket_2.tz
        @assert j1z + j3z == j2z + j4z "M quantum number not conserved in Hcm 2-body matrix element $(j1z) + $(j3z) != $(j2z) + $(j4z)"
        phase_bra = phase_ket = 1.0

        ## calculate phase factor due to fermion anticommutation
        phase_bra, phase_ket =  get_anticomm_phase(channel, 
                                                   bra1_bit, bra2_bit, ket1_bit, ket2_bit,
                                                   idx_m_bra_1, idx_m_bra_2, idx_m_ket_1, idx_m_ket_2)
        phase_factor = phase_bra * phase_ket
        Jmin_bra = div(abs(j1-j3), 2)
        Jmax_bra = div(abs(j1+j3), 2)
        Jmin_ket = div(abs(j2-j4), 2)
        Jmax_ket = div(abs(j2+j4), 2)

        Jmin = max(Jmin_bra, Jmin_ket)
        Jmax = min(Jmax_bra, Jmax_ket)

        key_1 = get_nkey4_shift(n1, l1, j1, t1)
        key_2 = get_nkey4_shift(n3, l3, j3, t3)
        key_3 = get_nkey4_shift(n2, l2, j2, t2)
        key_4 = get_nkey4_shift(n4, l4, j4, t4)

        for Jtot in Jmin:Jmax
            if abs(j1z + j3z) > Jtot * 2
                continue
            end
            nkey = (key_1, key_2, key_3, key_4, UInt(Jtot))
            tmp = dict_HCM[nkey]
            cg1 = dict_CGs[get_nkey6_shift(j1, j1z, j3, j3z, Jtot, div(j1z+j3z, 2); int_shift=int_shift)]
            cg2 = dict_CGs[get_nkey6_shift(j2, j2z, j4, j4z, Jtot, div(j2z+j4z, 2); int_shift=int_shift)]
            ret += tmp * cg1 * cg2 * phase_factor
        end
    end
    return ret
end

function delta_morb_except_m(ms1::SingleParticleState_Mscheme, ms2::SingleParticleState_Mscheme)
    if ms1.n == ms2.n && ms1.l == ms2.l && ms1.j == ms2.j && ms1.tz == ms2.tz
        return 1
    end
    return 0
end

"""
-◯●●◯- =(J_+)> -◯●◯●- =(J_-)> -◯●●◯-
"""
function eval_JmJp_diag(p_bitint_bra, n_bitint_bra, p_msps, n_msps)
    JmJpval = 0.0
    # not implemented yet
    for ch in ["p", "n"]
        msps = ifelse(ch == "p", p_msps, n_msps)
        bitint = ifelse(ch == "p", p_bitint_bra, n_bitint_bra)
        for idx_m in 1:length(msps) - 1
            onehot = Int128(1) << (idx_m-1)
            if (bitint & onehot) != onehot # if the bit is empty, skip
                continue
            end
            morb = msps[idx_m]
            # we are assuming msps are ordered by m, so we can just check the next one whether it is connectable
            idx_m_p1 = idx_m + 1
            next_onehot = Int128(1) << (idx_m_p1-1)
            if (bitint & next_onehot) == next_onehot # if the next bit is occupied, J+ cannot be applied
                continue
            end
            if delta_morb_except_m(morb, msps[idx_m_p1]) == 0
                continue
            end
            if msps[idx_m_p1].jz - morb.jz != 2
                continue # not connectable via J_+
            end
            jfac  = morb.j/2 * (morb.j/2 + 1) - msps[idx_m_p1].jz * morb.jz / 4
            #println("diag ch=$ch bra=ket=$(int2bitstr(bitint, length(msps))) => jfac = $jfac")
            JmJpval += jfac
        end
    end
    return JmJpval
end

"""
For hamming distance 4 (2 protons and 2 neutrons) cases, identifying the place of cre/ani operations are rather trivial.
However, for T=1 channels, one needs to be careful about 
"""
function eval_JmJp_nd(p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket, p_msps, n_msps)
    JmJpval = 0.0
    ## proton-neutron case
    if p_bitint_bra != p_bitint_ket && n_bitint_bra != n_bitint_ket
        p_cre = (p_bitint_bra ⊻ p_bitint_ket) & p_bitint_bra; idx_p_cre = max(0, log2(p_cre) .+ 1); idx_p_cre = Int(idx_p_cre)
        n_cre = (n_bitint_bra ⊻ n_bitint_ket) & n_bitint_bra; idx_n_cre = max(0, log2(n_cre) .+ 1); idx_n_cre = Int(idx_n_cre)
        p_ani = (p_bitint_bra ⊻ p_bitint_ket) & p_bitint_ket; idx_p_ani = max(0, log2(p_ani) .+ 1); idx_p_ani = Int(idx_p_ani)
        n_ani = (n_bitint_bra ⊻ n_bitint_ket) & n_bitint_ket; idx_n_ani = max(0, log2(n_ani) .+ 1); idx_n_ani = Int(idx_n_ani)
        j_p_cre = p_msps[idx_p_cre].j; j_n_cre = n_msps[idx_n_cre].j; m_p_cre = p_msps[idx_p_cre].jz; m_n_cre = n_msps[idx_n_cre].jz
        j_p_ani = p_msps[idx_p_ani].j; j_n_ani = n_msps[idx_n_ani].j; m_p_ani = p_msps[idx_p_ani].jz; m_n_ani = n_msps[idx_n_ani].jz
        if delta_morb_except_m(p_msps[idx_p_cre], p_msps[idx_p_ani]) == 0; return 0.0; end
        if delta_morb_except_m(n_msps[idx_n_cre], n_msps[idx_n_ani]) == 0; return 0.0; end
        if abs(m_p_cre - m_p_ani) != 2 || abs(m_n_cre - m_n_ani) != 2; return 0.0; end
        @assert m_p_cre - m_p_ani + m_n_cre - m_n_ani == 0 "Total angular momentum projection should be conserved, but got $(m_p_cre - m_p_ani) + $(m_n_cre - m_n_ani)"
        #println("proton: a^† = $(idx_p_cre) a = $(idx_p_ani)  neutron: a^† = $(idx_n_cre) a = $(idx_n_ani)")
        jfac = sqrt( j_p_cre/2 * (j_p_cre/2 + 1) - m_p_ani * m_p_cre/4 )
        jfac *= sqrt( j_n_cre/2 * (j_n_cre/2 + 1) - m_n_ani * m_n_cre/4 )
        JmJpval += jfac
    else # T=1 case
        ch = ifelse(p_bitint_bra != p_bitint_ket, "p", "n")
        msps = ifelse(ch == "p", p_msps, n_msps)
        bitint_bra = ifelse(ch == "p", p_bitint_bra, n_bitint_bra)
        bitint_ket = ifelse(ch == "p", p_bitint_ket, n_bitint_ket)
        cre_bits = (bitint_bra ⊻ bitint_ket) & bitint_bra
        ani_bits = (bitint_bra ⊻ bitint_ket) & bitint_ket

        # cre_bits and ani_bits should have 2 bits set, then, we need to get indices of the bits
        @assert count_ones(cre_bits) == 2 "There should be exactly 2 bits set in the creation bits, but got $(count_ones(cre_bits))"
        @assert count_ones(ani_bits) == 2 "There should be exactly 2 bits set in the annihilation bits, but got $(count_ones(ani_bits))"

        idx_cre_1 = idx_cre_2 = 0
        idx_ani_1 = idx_ani_2 = 0
        for idx_m in 1:length(msps)
            if (cre_bits & (Int128(1) << (idx_m-1))) != 0
                if idx_cre_1 == 0
                    idx_cre_1 = idx_m
                else
                    idx_cre_2 = idx_m
                end
            end
            if (ani_bits & (Int128(1) << (idx_m-1))) != 0
                if idx_ani_1 == 0
                    idx_ani_1 = idx_m
                else
                    idx_ani_2 = idx_m
                end
            end
        end

        # First, we assume idx_cre_1 - idx_ani_1, idx_cre_2 - idx_ani_2 are connectable with J^2
        if delta_morb_except_m(msps[idx_cre_1], msps[idx_ani_1]) == 0 
            idx_cre_1, idx_cre_2 = idx_cre_2, idx_cre_1 # swap the indices
        end

        if delta_morb_except_m(msps[idx_cre_2], msps[idx_ani_2]) == 0 
            return 0.0 # no longer connectable
        end

        j_cre_1 = msps[idx_cre_1].j; j_cre_2 = msps[idx_cre_2].j; m_cre_1 = msps[idx_cre_1].jz; m_cre_2 = msps[idx_cre_2].jz
        j_ani_1 = msps[idx_ani_1].j; j_ani_2 = msps[idx_ani_2].j; m_ani_1 = msps[idx_ani_1].jz; m_ani_2 = msps[idx_ani_2].jz
        if j_cre_1 != j_ani_1 || j_cre_2 != j_ani_2 || abs(m_cre_1 - m_ani_1) != 2 || abs(m_cre_2 - m_ani_2) != 2
            return 0.0
        end

        jfac  = sqrt( j_cre_1/2 * (j_cre_1/2 + 1) - m_ani_1 * m_cre_1/4 ) 
        jfac *= sqrt( j_cre_2/2 * (j_cre_2/2 + 1) - m_ani_2 * m_cre_2/4 ) 

        #println("nd T=1 channel $ch: bra $(int2bitstr(bitint_bra, length(msps))) ket $(int2bitstr(bitint_ket, length(msps)))  a^† = $(idx_cre_1), $(idx_cre_2) a = $(idx_ani_1), $(idx_ani_2) jfac $jfac")
        JmJpval += jfac
    end
    return JmJpval
end


"""
function construct_Jmat

Function to construct the J^2 operator, more specifically J_-J_+ to evaluate <J^2> = <J_-J_+> + Jz(Jz+1)
Since all the configurations for protons and neutrons are associated with the total angular momentum projection `Mp, Mn` and 
the `Mp+Mn` is to be conserved, one only needs to consider the first term of the J^2 operator, which is given via:

```math
\\hat{J}_{\\pm} = \\sum_{km} \\sqrt{j(j+1) - m(m\\pm 1)} a^{\\dagger}_{k,m\\pm 1} a_{k,m}
```
where `k` is the index of the single-particle state, `j` is the corresponding angular momentum, and `m` is the z-component of the angular momentum.

"""
function construct_Jmat(mdim, all_bitint_prod, p_msps, n_msps, verbose, to)
    Jmat = Dict{UInt, Float64}( )
    partial_dict = [ Dict{UInt, Float64}( ) for _ in 1:Threads.maxthreadid() ]
    @timeit to "thread loop" Threads.@threads for iter in 1:mdim^2
        idx_bra = div(iter-1, mdim) + 1
        idx_ket = iter - (idx_bra-1) * mdim
        if idx_bra > idx_ket
            continue 
        end

        p_bitint_bra, n_bitint_bra = all_bitint_prod[idx_bra]
        p_bitint_ket, n_bitint_ket = all_bitint_prod[idx_ket]

        p_ham_dist = count_ones( p_bitint_bra ⊻ p_bitint_ket )
        n_ham_dist = count_ones( n_bitint_bra ⊻ n_bitint_ket )

        sum_ham_dist = p_ham_dist + n_ham_dist
        if sum_ham_dist != 4 && sum_ham_dist != 0 # bra and ket cannot be connected by J^2 (still need to be checked!!!!!!)
            continue 
        end       
        
        tid = Threads.threadid()
        # # Calculate the <J_-J_+> matrix element
        Jij = 0.0
        if sum_ham_dist == 0
            Jij += eval_JmJp_diag(p_bitint_bra, n_bitint_bra, p_msps, n_msps)
        else
            Jij += eval_JmJp_nd(p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket, p_msps, n_msps)
        end
        if Jij == 0.0
            continue # no contribution to the J^2 operator
        end
        partial_dict[tid][hash_2ints(idx_bra, idx_ket)] = Jij
    end    

    # Combine the partial results from all threads
    @timeit to "collect" for tid in 1:Threads.maxthreadid()
        for (nkey, val) in partial_dict[tid]
            Jmat[nkey] = val
        end
    end

    if verbose >= 1 
        println("J matrix")
        for (nkey, val) in Jmat
            idx_bra, idx_ket = unhash_2ints(nkey)
            bra_bitint, ket_bitint = all_bitint_prod[idx_bra], all_bitint_prod[idx_ket]
            p_bra, n_bra = bra_bitint
            p_ket, n_ket = ket_bitint
            bra_bitstr = int2bitstr(p_bra, length(p_msps)) * " ⊗ " * int2bitstr(n_bra, length(n_msps))
            ket_bitstr = int2bitstr(p_ket, length(p_msps)) * " ⊗ " * int2bitstr(n_ket, length(n_msps))
            println("J[$(idx_bra),$(idx_ket)] = $(Jmat[nkey]) for bra $(bra_bitstr) ket $(ket_bitstr)")
         end
    end
    return Jmat
end

function eval_JJexpec(evecs, mdim, Jmat, Mtot, n_eigen, to)
    JJvals = zeros(Float64, n_eigen)
    Jzterm = Mtot/2 * (Mtot/2 + 1)
    v = zeros(Float64, mdim)
    Jv = zeros(Float64, mdim)
    partials = zeros(Float64, mdim, Threads.maxthreadid())
    keysvec = collect(keys(Jmat))
    for i in 1:n_eigen
        v .= evecs[:, i]
        operate_H_on_vec!(Jv, Jmat, v, partials, keysvec, to)
        JJvals[i] = dot(v, Jv) + Jzterm
    end
    return JJvals
end

function get_Jvals_fromJJ(JJvals; tol=1e-6) # J(J+1) = <J^2>. Note that J is doubled
    Jvals = zeros(Float64, length(JJvals))
    Jvals .= -1.0 # default value indicating not found
    for i in 1:length(JJvals)
        JJ = JJvals[i]
        for tJ2 in 0:1000
            Jtest = tJ2 / 2
            JJtest = Jtest * (Jtest + 1)
            if abs(JJtest - JJ) < tol
                Jvals[i] = Jtest
                break
            end
        end
    end
    return Jvals
end

function eval_H2(evecs, evals, Mask, FullHmat, org_mdim, to)
    evars = zeros(Float64, size(evecs, 2))
    psi = zeros(Float64, org_mdim)
    Hpsi = zeros(Float64, org_mdim)
    H2psi = zeros(Float64, org_mdim)
    partials = zeros(Float64, org_mdim, Threads.maxthreadid())
    keysvec = collect(keys(FullHmat))
    for i in 1:size(evecs, 2)
        evec = @view evecs[:, i]
        for key in keys(Mask)
            idx_full, idx_sub = unhash_2ints(key)
            psi[idx_full] = evec[idx_sub]
        end
        operate_H_on_vec!(Hpsi, FullHmat, psi, partials, keysvec, to)
        operate_H_on_vec!(H2psi, FullHmat, Hpsi, partials, keysvec, to)
        evars[i] = dot(psi, H2psi)
    end
    # If you only consider Hnd part, you don't need to subtract E^2
    return evars #.- evals.^2 
end