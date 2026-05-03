function _jj_to_ls_coeffs(orb_a::SingleParticleState, orb_b::SingleParticleState, J::Int; verbose::Bool=false)
    l_a = orb_a.l
    l_b = orb_b.l
    j2_a = orb_a.j
    j2_b = orb_b.j

    coeffs = Tuple{Int, Int, Float64}[]
    jhat_fac = sqrt((j2_a + 1) * (j2_b + 1))
    s = 1//2
    for S in 0:1
        for L in abs(l_a - l_b):(l_a + l_b)
            t9j = wigner9j(l_a, s, Rational(j2_a, 2),
                           l_b, s, Rational(j2_b, 2),
                           L,   S, J   )
            coef = jhat_fac * sqrt((2L + 1) * (2S + 1)) * t9j
            c = Float64(coef)
            if verbose
                println(" (L=$L, S=$S) with coeff=$(c)")
            end
            if abs(c) < 1.0e-14
                continue
            end
            push!(coeffs, (L, S, c))
        end
    end
    return coeffs
end


"""
Function to perform jj -> LS transformation using LS single-particle labels built by `define_ls_sps`.
For pn channel, key convention follows the input: (p,r)|(q,s).
"""
function _tbme_jj_to_ls_channel_by_lslabels(target, sps, ls_sps, dict_jj_to_ls, ch::Symbol; verbose::Bool=false)
    Vls = Dict{Tuple{UInt, UInt, Int, Int, Int, Int, Int}, Float64}()
    for (hash_i, hash_j) in keys(target)
        if ch == :pn
            p_jj, r_jj = unhash_2ints(hash_i)
            q_jj, s_jj = unhash_2ints(hash_j)

            p_ls = dict_jj_to_ls[p_jj]
            q_ls = dict_jj_to_ls[q_jj]
            r_ls = dict_jj_to_ls[r_jj]
            s_ls = dict_jj_to_ls[s_jj]
            hash_i_ls = hash_2ints(p_ls, r_ls)
            hash_j_ls = hash_2ints(q_ls, s_ls)

            orb_bra_1 = sps[p_jj]
            orb_bra_2 = sps[q_jj]
            orb_ket_1 = sps[r_jj]
            orb_ket_2 = sps[s_jj]

            for (Jf, vjj) in target[(hash_i, hash_j)]
                J = Int(round(Jf))
                bra_coeffs = _jj_to_ls_coeffs(orb_bra_1, orb_bra_2, J)
                ket_coeffs = _jj_to_ls_coeffs(orb_ket_1, orb_ket_2, J)
                for (Lbra, Sbra, cbra) in bra_coeffs
                    for (Lket, Sket, cket) in ket_coeffs
                        nkey = (hash_i_ls, hash_j_ls, J, Lbra, Sbra, Lket, Sket)
                        Vls[nkey] = get(Vls, nkey, 0.0) + cbra * cket * vjj
                    end
                end
            end
        else
            a_jj, b_jj = unhash_2ints(hash_i)
            c_jj, d_jj = unhash_2ints(hash_j)

            a_ls = dict_jj_to_ls[a_jj]
            b_ls = dict_jj_to_ls[b_jj]
            c_ls = dict_jj_to_ls[c_jj]
            d_ls = dict_jj_to_ls[d_jj]
            hash_i_ls = hash_2ints(a_ls, b_ls)
            hash_j_ls = hash_2ints(c_ls, d_ls)

            orb_a = sps[a_jj]
            orb_b = sps[b_jj]
            orb_c = sps[c_jj]
            orb_d = sps[d_jj]

            for (Jf, vjj) in target[(hash_i, hash_j)]
                J = Int(round(Jf))
                bra_coeffs = _jj_to_ls_coeffs(orb_a, orb_b, J)
                ket_coeffs = _jj_to_ls_coeffs(orb_c, orb_d, J)
                for (Lbra, Sbra, cbra) in bra_coeffs
                    for (Lket, Sket, cket) in ket_coeffs
                        nkey = (hash_i_ls, hash_j_ls, J, Lbra, Sbra, Lket, Sket)
                        Vls[nkey] = get(Vls, nkey, 0.0) + cbra * cket * vjj
                    end
                end
            end
        end
    end
    return Vls
end


function trans_orb_ab_to_label(orb_a::SingleParticleState, orb_b::SingleParticleState)
    ctz = ifelse(orb_a.tz == -1, "π", "ν")
    cl = chara_l[orb_a.l + 1]
    label = "$(ctz)$(orb_a.n)$(cl)$(orb_a.j)/2 "
    ctz = ifelse(orb_b.tz == -1, "π", "ν")
    cl = chara_l[orb_b.l + 1]
    label *= "$(ctz)$(orb_b.n)$(cl)$(orb_b.j)/2"
    return label
end


function eval_monopoleV(hamil::Hamiltonian_snt_fmt)
    sps = vcat(hamil.p_sps, hamil.n_sps)
    for ch in [:pp, :nn, :pn]
        target = nothing
        if ch == :pp
            target = hamil.V2b_pp
        elseif ch == :nn
            target = hamil.V2b_nn
        elseif ch == :pn
            target = hamil.V2b_pn
        end
        println("channel: $ch (#$(length(keys(target))))")
        dict_monopole = Dict{String, Float64}()
        for (hash_i, hash_j) in keys(target)
            if hash_i != hash_j && ch != :pn
                continue
            end            
            a, b = unhash_2ints(hash_i)
            if ch == :pn
                a, c = unhash_2ints(hash_i)
                b, d = unhash_2ints(hash_j)
                if a != c || b != d
                    continue
                end
            end
            if a > b
                continue
            end
            orb_a = sps[a]
            orb_b = sps[b]
            label = trans_orb_ab_to_label(orb_a, orb_b)
            vmono = 0.0
            Jsum = 0.0
            for (J, val) in target[hash_i, hash_j]
                vmono += (2J + 1) * val
                Jsum += 2J + 1
            end
            vmono /= Jsum
            dict_monopole[label] = vmono
            println("$label:", @sprintf("%12.6f", vmono))
        end    
    end
    return nothing
end


struct SingleParticleState_LS
    n::Int64
    l::Int64
    tz::Int64
end


function define_ls_sps(sps::Vector{SingleParticleState}, verbose::Bool=false)
    ls_sps = SingleParticleState_LS[ ]
    dict_jj_to_ls = Dict{Int, Int}()
    idx_ls = 1
    for idx_jj in 1:length(sps)
        sp = sps[idx_jj]
        n = sp.n
        l = sp.l
        j = sp.j                
        tz = sp.tz
        dict_jj_to_ls[idx_jj] = idx_ls
        if j == 2*l + 1
            # j = l + 1/2 -> assign to idx_ls
            idx_ls += 1
            push!(ls_sps, SingleParticleState_LS(n, l, tz))
        end
    end
    if verbose
        println("Defined LS-coupled single-particle states:")
        for idx_jj in 1:length(sps)
            sp = sps[idx_jj]
            n = sp.n
            l = sp.l
            j = sp.j                
            tz = sp.tz
            idx_ls = dict_jj_to_ls[idx_jj]
            ls_sp = ls_sps[idx_ls]
            println("jj idx: $idx_jj -> (n=$n, l=$l, j=$j, tz=$tz) -> ls idx: $idx_ls -> (n=$(ls_sp.n), l=$(ls_sp.l), tz=$(ls_sp.tz))")
        end
    end
    return ls_sps, dict_jj_to_ls
end


"""
Compute monopole values from jj-coupled spin-tensor decomposed channels.

For each two-body pair, this function performs a J-weighted average using
weight (2J+1) for each component of
(Total, central, LS, ALS, tensor).

Returns:
Dict(:pp => Dict{String, NTuple{5, Float64}},
     :nn => Dict{String, NTuple{5, Float64}},
     :pn => Dict{String, NTuple{5, Float64}})
"""
function eval_monopoleV_std(Vstd_channels::Dict, hamil::Hamiltonian_snt_fmt, ch::Symbol; verbose::Bool=true)
    sps = vcat(hamil.p_sps, hamil.n_sps)
    ret = Dict{Symbol, Dict{String, NTuple{5, Float64}}}()

    target = Vstd_channels
    num = Dict{Tuple{UInt, UInt}, NTuple{5, Float64}}()
    den = Dict{Tuple{UInt, UInt}, Float64}()
    label_of = Dict{Tuple{UInt, UInt}, String}()

    for ((hash_i, hash_j, J), vals) in target
        if ch != :pn && hash_i != hash_j
            continue
        end

        if ch == :pn
            a, c = unhash_2ints(hash_i)
            b, d = unhash_2ints(hash_j)
            if a != c || b != d
                continue
            end
            label_of[(hash_i, hash_j)] = trans_orb_ab_to_label(sps[a], sps[b])
        else
            a, b = unhash_2ints(hash_i)
            c, d = a, b
            label_of[(hash_i, hash_j)] = trans_orb_ab_to_label(sps[a], sps[b])
        end

        w = Float64(2 * J + 1)
        pkey = (hash_i, hash_j)
        old = get(num, pkey, (0.0, 0.0, 0.0, 0.0, 0.0))
        if abs(vals[1]) < 1.0e-14
            continue
        end
        num[pkey] = (
            old[1] + w * vals[1],
            old[2] + w * vals[2],
            old[3] + w * vals[3],
            old[4] + w * vals[4],
            old[5] + w * vals[5],
        )
        den[pkey] = get(den, pkey, 0.0) + w
        txt = "Total does not match sum of components for key ($hash_i, $hash_j, $J)\n"
        txt *= "<$a $b|V|$c $d>_J = $J -> vals=$(vals), w=$(w), num=$(num[pkey]), den=$(den[pkey])"    
        @assert num[pkey][1] ≈ num[pkey][2] + num[pkey][3] + num[pkey][4] + num[pkey][5] "$txt"
    end

    mono_ch = Dict{String, NTuple{5, Float64}}()
    for (pkey, nvals) in num
        wsum = get(den, pkey, 0.0)
        if wsum <= 0.0
            continue
        end
        label = label_of[pkey]
        mono_ch[label] = (
            nvals[1] / wsum,
            nvals[2] / wsum,
            nvals[3] / wsum,
            nvals[4] / wsum,
            nvals[5] / wsum,
        )
        if verbose
            v = mono_ch[label]
            println("channel $ch | $label | ",
                    @sprintf("%12.8f", v[1]),
                    @sprintf("%12.8f", v[2]),
                    @sprintf("%12.8f", v[3]),
                    @sprintf("%12.8f", v[4]),
                    @sprintf("%12.8f", v[5]))
        end
    end
    ret[ch] = mono_ch

    return ret
end


function _std_Kirson_channel(target, sps, ls_sps, dict_jj_to_ls, ch::Symbol;
                             show_breakdown::Bool=false)
    Vls_label = _tbme_jj_to_ls_channel_by_lslabels(target, sps, ls_sps, dict_jj_to_ls, ch)
    Vstd = Dict{Tuple{UInt, UInt, Int}, NTuple{5, Float64}}()

    for ((hash_i, hash_j), vecJ) in target
        if ch == :pn
            p, r = unhash_2ints(hash_i)
            q, s = unhash_2ints(hash_j)
            orb_bra_1 = sps[p]
            orb_bra_2 = sps[q]
            orb_ket_1 = sps[r]
            orb_ket_2 = sps[s]
            hash_i_ls = hash_2ints(dict_jj_to_ls[p], dict_jj_to_ls[r])
            hash_j_ls = hash_2ints(dict_jj_to_ls[q], dict_jj_to_ls[s])
        else
            p, q = unhash_2ints(hash_i)
            r, s = unhash_2ints(hash_j)
            orb_bra_1 = sps[p]
            orb_bra_2 = sps[q]
            orb_ket_1 = sps[r]
            orb_ket_2 = sps[s]
            hash_i_ls = hash_2ints(dict_jj_to_ls[p], dict_jj_to_ls[q])
            hash_j_ls = hash_2ints(dict_jj_to_ls[r], dict_jj_to_ls[s])
        end

        for (Jf, _) in vecJ
            J = Int(round(Jf))
            bra_coeffs = _jj_to_ls_coeffs(orb_bra_1, orb_bra_2, J)
            ket_coeffs = _jj_to_ls_coeffs(orb_ket_1, orb_ket_2, J)

            central = 0.0
            ls_term = 0.0
            als_term = 0.0
            tensor = 0.0

            for (Lbra, Sbra, cbraJ) in bra_coeffs
                for (Lket, Sket, cketJ) in ket_coeffs
                    for rank in 0:2
                        s6_out = wigner6j(Float64, Lbra, Sbra, J, Sket, Lket, rank)
                        abs(s6_out) < 1.0e-14 && continue

                        Jp_min = max(abs(Lbra - Sbra), abs(Lket - Sket))
                        Jp_max = min(Lbra + Sbra, Lket + Sket)
                        if Jp_min > Jp_max
                            continue
                        end

                        inner = 0.0
                        for Jp in Jp_min:Jp_max
                            nkey = (hash_i_ls, hash_j_ls, Jp, Lbra, Sbra, Lket, Sket)
                            vls_jp = get(Vls_label, nkey, 0.0)
                            abs(vls_jp) < 1.0e-14 && continue
                            s6_in = wigner6j(Float64, Lbra, Sbra, Jp, Sket, Lket, rank)
                            abs(s6_in) < 1.0e-14 && continue
                            inner += ((-1.0)^Jp) * Float64(2 * Jp + 1) * s6_in * vls_jp
                        end
                        abs(inner) < 1.0e-14 && continue

                        contrib = ((-1.0)^J) * Float64(2 * rank + 1) * cbraJ * cketJ * s6_out * inner
                        if rank == 0
                            central += contrib
                        elseif rank == 1
                            if Sbra == Sket
                                ls_term += contrib
                            else
                                als_term += contrib
                            end
                        else
                            tensor += contrib
                        end
                    end
                end
            end

            total = central + ls_term + als_term + tensor
            Vstd[(hash_i, hash_j, J)] = (total, central, ls_term, als_term, tensor)
            canonically_orderd = p < r || (p == r && q <= s)
            if show_breakdown && canonically_orderd
                println("$ch <$p $q|V|$r $s>_J = $J -> ",
                        @sprintf("%12.6f", total),
                        @sprintf("%12.6f", central),
                        @sprintf("%12.6f", ls_term),
                        @sprintf("%12.6f", als_term),
                        @sprintf("%12.6f", tensor))
            end
        end
    end

    # Evaluating monopole values from the spin-tensor decomposed channels
    # This is achieved by taking J-weighted average of each component, with weight (2J+1) for each J value.
    num = Dict{Tuple{UInt, UInt}, NTuple{5, Float64}}()
    den = Dict{Tuple{UInt, UInt}, Float64}()
    label_of = Dict{Tuple{UInt, UInt}, String}()

    for ((hash_i, hash_j, J), vals) in Vstd
        if ch != :pn && hash_i != hash_j
            continue
        end
        if ch == :pn
            a, c = unhash_2ints(hash_i)
            b, d = unhash_2ints(hash_j)
            if a != c || b != d
                continue
            end
            label_of[(hash_i, hash_j)] = trans_orb_ab_to_label(sps[a], sps[b])
        else
            a, b = unhash_2ints(hash_i)
            c, d = a, b
            label_of[(hash_i, hash_j)] = trans_orb_ab_to_label(sps[a], sps[b])
        end
        if a > b
            continue
        end

        w = Float64(2 * J + 1)
        pkey = (hash_i, hash_j)
        old = get(num, pkey, (0.0, 0.0, 0.0, 0.0, 0.0))
        if abs(vals[1]) < 1.0e-14
            continue
        end
        num[pkey] = (
            old[1] + w * vals[1],
            old[2] + w * vals[2],
            old[3] + w * vals[3],
            old[4] + w * vals[4],
            old[5] + w * vals[5],
        )
        den[pkey] = get(den, pkey, 0.0) + w
        txt = "Total does not match sum of components for key ($hash_i, $hash_j, $J)\n"
        txt *= "<$a $b|V|$c $d>_J = $J -> vals=$(vals), w=$(w), num=$(num[pkey]), den=$(den[pkey])"    
        @assert num[pkey][1] ≈ num[pkey][2] + num[pkey][3] + num[pkey][4] + num[pkey][5] "$txt"
    end

    mono_ch = Dict{String, NTuple{5, Float64}}()
    for (pkey, nvals) in num
        wsum = get(den, pkey, 0.0)
        if wsum <= 0.0
            continue
        end
        label = label_of[pkey]
        mono_ch[label] = (
            nvals[1] / wsum,
            nvals[2] / wsum,
            nvals[3] / wsum,
            nvals[4] / wsum,
            nvals[5] / wsum,
        )
        v = mono_ch[label]
    
        println("channel $ch | $label: | ",
                @sprintf("%12.8f", v[1]),
                @sprintf("%12.8f", v[2]),
                @sprintf("%12.8f", v[3]),
                @sprintf("%12.8f", v[4]),
                @sprintf("%12.8f", v[5]))
    end
   

    return Vstd
end

"""
Compute directly the spin-tensor decomposition in jj-coupling with Kirson's formulation.
The essential thing is to sum over all possible a', b', c', d' and J' intermediate states
mediated by the relation through L, S values.

<ab; JT|V|cd; JT> = (-1)^J  (2k+1) \\sum_{L,Lp,S,Sp} (ab|LSJ) (cd|L'S'J) {L S J; S' L' k} \\sum_{Jp} (-1)^{Jp} (2Jp+1)
× {L S Jp; S' L' k} \\sum_{ja', jb', jc', jd'} (ja' jb'|LSJp) (jc' jd'|L'S'Jp) <ja' jb'; Jp|V|jc' jd'; Jp>
"""
function std_Kirson(hamil::Hamiltonian_snt_fmt; verbose::Bool=false)
    sps = vcat(hamil.p_sps, hamil.n_sps)
    ls_sps, dict_jj_to_ls = define_ls_sps(sps, false)
    Vstd_pp = _std_Kirson_channel(hamil.V2b_pp, sps, ls_sps, dict_jj_to_ls, :pp; show_breakdown=verbose)
    Vstd_nn = _std_Kirson_channel(hamil.V2b_nn, sps, ls_sps, dict_jj_to_ls, :nn; show_breakdown=verbose)
    Vstd_pn = _std_Kirson_channel(hamil.V2b_pn, sps, ls_sps, dict_jj_to_ls, :pn; show_breakdown=verbose)
    return Dict(:pp => Vstd_pp, :nn => Vstd_nn, :pn => Vstd_pn)
end


"""
spin_tensor_decomposition(hamil::Hamiltonian_snt_fmt)

Function to decompose the two-body matrix elements into the rank of sphericakl tensor.
This is achieved by firstly tranforming the TBMEs from jj-coupling to LS-coupling,
and then performing the decomposition.
"""
function spin_tensor_decomposition(hamil::Hamiltonian_snt_fmt;
                                   verbose::Bool=false,
                                   is_check_roundtrip::Bool=true)

    sps = vcat(hamil.p_sps, hamil.n_sps)
    LS_sps, dict_jj_to_ls = define_ls_sps(sps)
    std_Kirson(hamil; verbose=verbose)
   
    return true
end