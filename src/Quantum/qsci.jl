struct Mconfig
    Nocc::Int
    Nocc_p::Int
    Nocc_n::Int
    Nexc::Int
    Mtot::Int
    Mp::Int
    Mn::Int
    parity::Int
end

"""
Read a p-n combined bitstring, and evaluate Nocc, parity, 
"""
function read_bitstr(bitstr::String, m_sps::Vector{SingleParticleState_Mscheme};
                     Nexc_ref=0, Qiskit_ordered::Bool=true) ::Mconfig
    bitstr = replace(bitstr, "⊗"=> "") 
    bitstr = replace(bitstr, " "=> "") 
    bitstr = ifelse(Qiskit_ordered, reverse(bitstr), bitstr)

    Nq = length(bitstr)
    @assert Nq == length(m_sps) "Length of bitstr must match the number of single-particle states"
    Nocc = Nocc_p = Nocc_n = Mtot = Mtot_p = Mtot_n = 0
    Nexc = -Nexc_ref
    parity = 1
    for i = 1:Nq
        if bitstr[i] == '1'
            orbit = m_sps[i]
            e = orbit.e
            l = orbit.l
            parity *=(-1)^l
            Nocc += 1
            if orbit.tz == 1
                Nocc_p += 1
            elseif orbit.tz == -1
                Nocc_n += 1
            else
                error("Invalid tz value: $(orbit.tz)")
            end
            Nexc += e
            Mtot += orbit.jz
            if orbit.tz == -1
                Mtot_p += orbit.jz
            elseif orbit.tz == 1
                Mtot_n += orbit.jz
            else
                error("Invalid tz value: $(orbit.tz)")
            end
        end
    end
    return Mconfig(Nocc, Nocc_p, Nocc_n, Nexc, Mtot, Mtot_p, Mtot_n, parity)
end

function read_bitint(bint, m_sps::Vector{SingleParticleState_Mscheme}; Nexc_ref=0) ::Mconfig
    parity = 1
    Nocc = Nocc_p = Nocc_n = Mtot = Mtot_p = Mtot_n = Nexc = 0
    for idx in 1:length(m_sps)
        if (bint & (Int128(1) << (idx-1))) != 0
            parity *= (-1)^(m_sps[idx].l)
            m = m_sps[idx].jz
            Mtot += m
            Nocc += 1
            Nexc += m_sps[idx].e
            if m_sps[idx].tz == -1
                Mtot_p += m
                Nocc_p += 1
            elseif m_sps[idx].tz == 1
                Mtot_n += m
                Nocc_n += 1
            end
        end
    end
    Nexc = Nexc - Nexc_ref
    return Mconfig(Nocc, Nocc_p, Nocc_n, Nexc, Mtot, Mtot_p, Mtot_n, parity)
end

struct Hamiltonian_snt_fmt
    lp::Int # number of proton single-particle states
    ln::Int # number of neutron single-particle states
    cp::Int # number of proton core 
    cn::Int # number of neutron core
    p_sps::Vector{SingleParticleState} # proton single-particle states
    n_sps::Vector{SingleParticleState} # neutron single-particle states
    h_1b::Dict{Tuple{Int,Int}, Float64} # one-body matrix elements
    V2b_pp::Dict{Tuple{UInt,UInt}, Vector{Tuple{Float64,Float64}}} # proton-proton two-body matrix elements; 1st UInt = (p, q), 2nd UInt = (r, s)
    V2b_nn::Dict{Tuple{UInt,UInt}, Vector{Tuple{Float64,Float64}}} # neutron-neutron two-body matrix elements; like V2b_pp
    V2b_pn::Dict{Tuple{UInt,UInt}, Vector{Tuple{Float64,Float64}}} # proton-neutron two-body matrix elements; 1st UInt = (p, r), 2nd UInt = (q, s)
    V3b::Vector{Obj_3BME}
end

function sps2msps(sps::Vector{SingleParticleState})
    m_sps = SingleParticleState_Mscheme[ ]
    for orb in sps
        n = orb.n
        l = orb.l
        j = orb.j
        tz = orb.tz
        e = 2n + l
        for jz in -j:2:j
            push!(m_sps, SingleParticleState_Mscheme(e, n, l, j, jz, tz))
        end
    end
    return m_sps
end

function make_m2j_j2m_dicts(p_sps::SPS, n_sps::SPS) where SPS <: Vector{SingleParticleState}
    dict_m2j = Dict{Int, Int}( )
    dict_j2m = Dict{Int, Vector{Int}}( )
    idx_m = 0
    for (idx_j, sps) in enumerate(p_sps)
        j = sps.j
        dict_j2m[idx_j] = Int[ ]
        for jz in -j:2:j
            idx_m += 1
            dict_m2j[idx_m] = idx_j
            push!(dict_j2m[idx_j], idx_m)
        end
    end
    for (idx_j, sps) in enumerate(n_sps)
        j = sps.j
        dict_j2m[idx_j + length(p_sps)] = Int[ ]
        for jz in -j:2:j
            idx_m += 1
            dict_m2j[idx_m] = idx_j + length(p_sps) # offset by the number of proton states
            push!(dict_j2m[idx_j + length(p_sps)], idx_m)
        end
    end
    return dict_m2j, dict_j2m
end

function hash_2ints(i::Int, j::Int) ::UInt
    return UInt(i) << 32 | UInt(j)
end

function hash_2ints(i::Int128, j::Int128) ::UInt128
    return UInt128(i) << 64 | UInt128(j)
end

function unhash_2ints(h::UInt) :: Tuple{Int, Int}
    i = Int(h >> 32)
    j = Int(h & 0xFFFFFFFF)
    return (i, j)
end

function read_smsnt_dev(sntf, 
                        cnuc::String;
                        ncsm::Bool=false, 
                        reorder_orbits::Bool=false,
                        ignore_massop=false) 
    idx_dicts_org2reo = Dict{Int, Int}()
    Anum = parse(Int64,match(reg,cnuc).match)
    if Anum <= 0
        @error "Invalid nucleus $cnuc, Anum=$Anum"
        exit()
    end
    f = open(sntf,"r");tlines = readlines(f);close(f)
    lines = rm_comment(tlines)
    line = lines[1]    
    lp,ln,cp,cn = map(x->parse(Int,x),rm_nan(split(line," ")))
    p_sps = SingleParticleState[ ]
    n_sps = SingleParticleState[ ]
    sps_both = SingleParticleState[ ]
    for i = 1:lp
        ith, n, l, j, tz = map(x->parse(Int64,x),rm_nan(split(lines[1+i]," "))[1:5])
        push!(p_sps, SingleParticleState(2*n+l, n, l, j, tz, [0.0], [false], [true], [false]))
        push!(sps_both, SingleParticleState(2*n+l, n, l, j, tz, [0.0], [false], [true], [false]))
    end
    for i = 1:ln
        ith, n, l, j, tz = map(x->parse(Int64,x),rm_nan(split(lines[1+i+ln]," "))[1:5])
        push!(n_sps, SingleParticleState(2*n+l, n, l, j, tz, [0.0], [false], [true], [false]))
        push!(sps_both, SingleParticleState(2*n+l, n, l, j, tz, [0.0], [false], [true], [false]))  
    end
    if reorder_orbits
        # Reorder the orbits by j (smaller j first), and then by n (larger n first) if j is the same. 
        # First we need to update the idx_dicts_org2reo to reflect the new ordering,
        # which will be used later when we read the 1b and 2b matrix elements.
        idxs_p = sortperm(p_sps, by = x -> (x.j, -x.n))
        idxs_n = sortperm(n_sps, by = x -> (x.j, -x.n))
        for (new_idx, old_idx) in enumerate(idxs_p)
            idx_dicts_org2reo[old_idx] = new_idx
        end
        for (new_idx, old_idx) in enumerate(idxs_n)
            idx_dicts_org2reo[old_idx + lp] = new_idx + lp
        end        
        p_sps = sort(p_sps, by = x -> (x.j, -x.n))
        n_sps = sort(n_sps, by = x -> (x.j, -x.n))
        sps_both = vcat(p_sps, n_sps)
    end

    # num_of_1b-term, method
    nsp, zero = map(x->parse(Int,x),rm_nan(split(lines[1+ln+lp+1]," "))[1:2])
    
    # one-body matrix elements; For ones within 1hw space, only diagonal (single-particle energies) components appear.
    # Something like: 
    # 1    1      0.75000000
    # 4    1      0.61237244
    # 2    2      1.25000000
    h_1b = Dict{Tuple{Int,Int}, Float64}()
    for n = 1:nsp
        ttxt = rm_nan(split(lines[1+ln+lp+1+n]," "))
        i, j = map(x->parse(Int,x),ttxt[1:2])        
        spe = parse(Float64,ttxt[3])  
        if reorder_orbits
            i = idx_dicts_org2reo[i]
            j = idx_dicts_org2reo[j]
        end
        if i == j
            h_1b[(i, j)] = spe 
        else
            h_1b[(i, j)] = spe 
            h_1b[(j, i)] = spe 
        end
    end
    # The next line would be...  6572 10 20.00 (NuHamil case) or 158   1  18 -0.300000 (KSHELL snt)
    # NuHamil: num_of_TBME, hw_option, hw-value
    # KSHELL:  num_of_TBME, mass_option, Aref, power for mass dependence
    tmp = rm_nan(split(lines[1+ln+lp+1+nsp+1]," "))
    @assert 3 <= length(tmp) <= 4 "Invalid snt file format: expected 3 or 4 elements in the line after 1b terms, got $(length(tmp)) elements. Check the file format."
    snt_mode = ifelse(length(tmp) == 3 && tmp[2]!="0", "NuHamil", "KSHELL")

    ntbme = hw_option = hw_value = massop = Aref = p = 0
    if snt_mode == "NuHamil"
        ntbme, hw_option, hw_value = tmp[1:3]
        ntbme = parse(Int,ntbme); hw_option = parse(Int,hw_option); hw_value = parse(Float64,hw_value)
        Aref = Anum
        println("ncsm $ncsm, ntbme $ntbme hw_option $hw_option, hw_value $hw_value")
        if !ncsm 
            @warn "ncsm is set false, but the interaction file is in NuHamil format. I would suggest to check it."
        end
    else
        p = 1.0
        try
            ntbme,massop,Aref,p = tmp[1:4]
            p = parse(Float64,p)
        catch
            ntbme, massop, Aref = tmp[1:3]
        end
        ntbme = parse(Int,ntbme); massop=parse(Int,massop)
        Aref = parse(Float64,string(Aref))
    end
    Adep_factor_1b = ifelse(ncsm, hw_value *(Anum-1)/Anum, 1.0)
    if ncsm && hw_option == 10 && !ignore_massop
        for i in 1:nsp
            for j in i:nsp
                if haskey(h_1b, (i, j))
                    jfac =sqrt( 1 + sps_both[i].j) # no need
                    h_1b[(i, j)] *= Adep_factor_1b
                    if i != j 
                        h_1b[(j, i)] *= Adep_factor_1b 
                    end
                end
            end
        end
    end
    vfactor = ifelse(massop==1, (Anum/Aref)^p, 1.0)

    # TBME lines start from now on
    # (K)    1   1   1   1    0    -1.6913000
    # (N)    1   1   1   1    0    -8.22398567     -0.00000000 (vpp term follows)
    # ---
    # Vpp (Vnn is similar)
    # V_{pqrs} π^†_p π^†_q π_s π_r 
    V2b_pp = Dict{Tuple{UInt, UInt}, Vector{Tuple{Float64, Float64}}}( )
    V2b_nn = Dict{Tuple{UInt, UInt}, Vector{Tuple{Float64, Float64}}}( )
    V2b_pn = Dict{Tuple{UInt, UInt}, Vector{Tuple{Float64, Float64}}}( )
    println("massop $massop hw_option $hw_option vfactor $vfactor Anum $Anum Aref $Aref")
    for ith = 1:ntbme
        tmp = rm_nan(split(lines[1+ln+lp+1+nsp+1+ith], " "))
        p, q, r, s, totJ, TBME = tmp[1:6]
        Vkin = 0.0
        if length(tmp) > 6
            try
                Vkin = parse(Float64, split(tmp[7]))  # vpp term, if present
            catch
                Vkin = 0.0
            end
        end
        p = parse(Int,p); q = parse(Int,q); r = parse(Int,r); s = parse(Int,s)
        if reorder_orbits
            p = idx_dicts_org2reo[p]
            q = idx_dicts_org2reo[q]
            r = idx_dicts_org2reo[r]
            s = idx_dicts_org2reo[s]
        end
        totJ = parse(Float64,totJ)
        TBME = parse(Float64,TBME)

        # Check the ordering of indices and apply the necessary phase factor to make them canonically ordered.
        j_p = sps_both[p].j; j_q = sps_both[q].j; j_r = sps_both[r].j; j_s = sps_both[s].j
        if p > q
            p, q = q, p
            TBME *= (-1)^(div(j_p+j_q, 2) + totJ + 1)
        end
        if r > s
            r, s = s, r
            TBME *= (-1)^(div(j_r+j_s, 2) + totJ + 1)
        end
        if p > r || (p == r && q > s)
            p, q, r, s = r, s, p, q
        end

        ## Mass dependence factor if needed
        if massop == 1 && !ignore_massop
            TBME *= vfactor
        end
        if hw_option == 10 && !ignore_massop
            TBME = TBME + Vkin * hw_value / Anum
        end
        pnrank = 0
        if p <= lp && q <= lp
            pnrank = 1 # proton-proton
            TotT = 1
        elseif p > lp && q > lp
            pnrank = 3 # neutron-neutron
            TotT = 1
        elseif p <= lp && q > lp
            pnrank = 2
            TotT = 0
        else
            @error "This should not happen in a valid snt file: p=$p, q=$q, lp=$lp, ln=$ln"
            exit()
        end
        TotT = abs(pnrank - 2)

        if TotT == 1
            bra_Uint = hash_2ints(p, q)
            ket_Uint = hash_2ints(r, s)
            phase_bra = phase_ket = 1.0
            if p > q
                bra_Uint = hash_2ints(q, p)
                phase_bra *= (-1)^(div(sps_both[p].j+sps_both[q].j, 2) + totJ + 1)
            end
            if r > s
                ket_Uint = hash_2ints(s, r)
                phase_ket *= (-1)^(div(sps_both[r].j+sps_both[s].j, 2) + totJ + 1)
            end
            target = ifelse(pnrank == 1, V2b_pp, V2b_nn)
            if !haskey(target, (bra_Uint, ket_Uint))
                target[(bra_Uint, ket_Uint)] = Vector{Tuple{Float64, Float64}}( )
            end
            push!(target[(bra_Uint, ket_Uint)], (totJ, TBME*phase_bra*phase_ket))
            #println("p $p q $q r $r s $s totJ $totJ TBME $TBME key ($bra_Uint, $ket_Uint) ")
            # permutation for bra and ket
            if bra_Uint != ket_Uint
                if !haskey(target, (ket_Uint, bra_Uint))
                    target[(ket_Uint, bra_Uint)] = Vector{Tuple{Float64, Float64}}( )
                end
                push!(target[(ket_Uint, bra_Uint)], (totJ, TBME*phase_bra*phase_ket))
            end

        elseif TotT == 0
            @assert p < q && r < s "Indices in bra/ket must be ordered for TotT=0 case: p=$p, q=$q, r=$r, s=$s"
            p_Uint = hash_2ints(p, r)
            n_Uint = hash_2ints(q, s)
            if !haskey(V2b_pn, (p_Uint, n_Uint))
                V2b_pn[(p_Uint, n_Uint)] = Vector{Tuple{Float64, Float64}}( )
            end
            push!(V2b_pn[(p_Uint, n_Uint)], (totJ, TBME))

            if p != r || q != s
                p_Uint = hash_2ints(r, p)
                n_Uint = hash_2ints(s, q)
                if !haskey(V2b_pn, (p_Uint, n_Uint))
                    V2b_pn[(p_Uint, n_Uint)] = Vector{Tuple{Float64, Float64}}( )
                end
                push!(V2b_pn[(p_Uint, n_Uint)], (totJ, TBME))
            end
        else
            @error "Invalid TotT value: $TotT"
        end
    end
    println("# of entries 1b terms: $(length(h_1b)) Vpp $(length(V2b_pp)) Vnn $(length(V2b_nn)) Vpn $(length(V2b_pn))")
    # Prepare undefined V3b
    V3b = Obj_3BME[]
    return Hamiltonian_snt_fmt(lp, ln, cp, cn, p_sps,n_sps, h_1b, V2b_pp, V2b_nn, V2b_pn, V3b)
end


"""

Evaluate the Hamiltonian matrix element H_ij for some non-diagonal cases.
The 02&20 means that either the proton or neutron part have hamming distance of 2,
and the other part has hamming distance of 0.
(0, 2) => n1b, nn, pn, nnn, ppn, pnn, (2, 0) => p1b, pp, pn, ppp, ppn, pnn

Note: The way to evaluate phase factor coming from anti-commutation relations here is a bit tricky.
Since we already know one of the components in both bra and ket is identical, which is referred to as 
spectator in the code.
Let us give an example: A configuration in 8He within p shell, < 000000 ⊗ 100111 | H | 000000 ⊗ 101110 >.
Under the operators with (0, 2) hamming distances, 2nd, 3rd, and 6th neutron orbits (counting from the right most bit),
can be a spectator.
The involved operators will be a^†_1 a^†_2 a_2 a_4, a^†_1 a^†_3 a_3 a_4, a^†_1 a^†_6 a_4 a_6, respectively.
We firstly evaluate the phase factor coming from the "target"
"""
function eval_Hij_nondiag_02_20(Hamil_snt, p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket,
                                p_ham_dist, n_ham_dist, p_msps, n_msps, dict_m2j, Hrank, int_shift, dict_CGs, verbose)
    H_ij = 0.0
    # proton/neutron 1-body term (irrelevant for model space smaller than two major shells)
    ch = ifelse(p_ham_dist == 2, "p", "n")

    m_sps = ifelse(ch == "p", p_msps, n_msps)
    bitint_bra = ifelse(ch == "p", p_bitint_bra, n_bitint_bra)
    bitint_ket = ifelse(ch == "p", p_bitint_ket, n_bitint_ket)
    offset_idx = ifelse(ch == "p", 0, length(p_msps))

    cre = ani = 0
    phase_bra = phase_ket = 1
    for idx_m in 1:length(m_sps)
        onehot = Int128(1) << (idx_m-1)
        if (bitint_ket & onehot) == onehot && (bitint_bra & onehot) == 0
            ani = idx_m
        end
        if (bitint_bra & onehot) == onehot && (bitint_ket & onehot) == 0
            cre = idx_m
        end
    end
    @assert ani != 0 && cre != 0 "ani or cre must be nonzero: ani = $ani, cre = $cre bitint_bra $bitint_bra bitint_ket $bitint_ket"

    for idx_m in 1:cre-1
        onehot = Int128(1) << (idx_m-1)
        if (bitint_ket & onehot) == onehot 
            phase_bra *= -1
        end
    end     
    for idx_m in 1:ani-1
        onehot = Int128(1) << (idx_m-1)
        if (bitint_bra & onehot) == onehot
            phase_ket *= -1
        end
    end 
    if verbose >= 1
        println("cre $cre ani $ani phase_bra $phase_bra phase_ket $phase_ket")
    end
    
    idx_jj_ani = dict_m2j[ani + offset_idx]
    idx_jj_cre = dict_m2j[cre + offset_idx]
    tkey = (idx_jj_cre, idx_jj_ani)
    if haskey(Hamil_snt.h_1b, tkey)        
        if verbose >=1
            println("H_ij $(Hamil_snt.h_1b[tkey]): 1b term from ($idx_jj_cre, $idx_jj_ani)")
        end
        H_ij += Hamil_snt.h_1b[tkey] * phase_bra * phase_ket  *(-1) 
    end
    # Two-body terms
    # Other than ani/cre bit, the rest of the occupied bits must be the same, which can be called "spectators"
    # One should run through all the spectators to sum up the two-body matrix elements
    if Hrank >= 2
        hamil = ifelse(ch == "p", Hamil_snt.V2b_pp, Hamil_snt.V2b_nn)
        for idx_spec in 1:length(m_sps)
            if idx_spec == ani || idx_spec == cre
                continue 
            end
            phase = phase_bra * phase_ket * (-1) # (-1) from to make spectator ops in vicinity

            idx_jj_spec = dict_m2j[idx_spec + offset_idx]
            onehot_spec = Int128(1) << (idx_spec-1)
            if (bitint_bra & onehot_spec) != onehot_spec || (bitint_ket & onehot_spec) != onehot_spec
                continue 
            end
            bra_p, bra_q = cre, idx_spec
            if cre > idx_spec
                bra_p, bra_q = idx_spec, cre
                phase *= -1
            end
            bra_p_jj = dict_m2j[bra_p + offset_idx]
            bra_q_jj = dict_m2j[bra_q + offset_idx]
            bra_Uint_jj = hash_2ints(bra_p_jj, bra_q_jj)

            ket_r, ket_s = ani, idx_spec
            if ani > idx_spec
                ket_r, ket_s = idx_spec, ani
                phase *= -1
            end            
            ket_r_jj = dict_m2j[ket_r + offset_idx]
            ket_s_jj = dict_m2j[ket_s + offset_idx]
            ket_Uint_jj = hash_2ints(ket_r_jj, ket_s_jj)
            if verbose >= 1
                println("< $bra_p $bra_q| $ket_r $ket_s >  =(jj)=> <$bra_p_jj $bra_q_jj| $ket_r_jj $ket_s_jj>")
            end
            if haskey(hamil, (bra_Uint_jj, ket_Uint_jj))
                Vjj_relevant = hamil[(bra_Uint_jj, ket_Uint_jj)]
                j1 = m_sps[bra_p].j; m1 = m_sps[bra_p].jz
                j2 = m_sps[bra_q].j; m2 = m_sps[bra_q].jz
                j3 = m_sps[ket_r].j; m3 = m_sps[ket_r].jz
                j4 = m_sps[ket_s].j; m4 = m_sps[ket_s].jz
                Mbra = div(m1+m2, 2) 
                Mket = div(m3+m4, 2)
                for (J, V_jj) in Vjj_relevant
                    J = Int(J)
                    if abs(Mbra) > J; continue; end
                    if abs(Mket) > J; continue; end
                    NJ_bra = sqrt(deltaf(idx_jj_cre, idx_jj_spec) *(-1)^J + 1.0) / (1+deltaf(idx_jj_cre, idx_jj_spec))
                    NJ_ket = sqrt(deltaf(idx_jj_ani, idx_jj_spec) *(-1)^J + 1.0) / (1+deltaf(idx_jj_ani, idx_jj_spec))
                    cg1 = dict_CGs[get_nkey6_shift(j1, m1, j2, m2, J, Mbra; int_shift=int_shift)]
                    cg2 = dict_CGs[get_nkey6_shift(j3, m3, j4, m4, J, Mket; int_shift=int_shift)]
                    CG = cg1 * cg2
                    vtmp = V_jj * CG / (NJ_bra*NJ_ket) *phase
                    if verbose >= 1
                        println("Vjj(02/20): J=$J v=$V_jj vtmp $vtmp CG $CG phase $phase")
                    end
                    H_ij += vtmp
                end
            end
        end
    end

    return H_ij 
end

"""
a^†_p a^†_q a_r a_s where {p, q} ^ {r, s} = ∅
"""
function eval_Hij_nondiag_04_40(Hamil_snt, p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket,
                                p_ham_dist, n_ham_dist, p_msps, n_msps, dict_m2j, Hrank, int_shift, dict_CGs, verbose)
    H_ij = 0.0
    ch = ifelse(p_ham_dist == 4, "p", "n")
    m_sps = ifelse(ch == "p", p_msps, n_msps)
    bitint_bra = ifelse(ch == "p", p_bitint_bra, n_bitint_bra)
    bitint_ket = ifelse(ch == "p", p_bitint_ket, n_bitint_ket)
    offset_idx = ifelse(ch == "p", 0, length(p_msps))
    hamil = ifelse(ch == "p", Hamil_snt.V2b_pp, Hamil_snt.V2b_nn)

    cre1 = cre2 = ani1 = ani2 = 0
    for idx_m in 1:length(m_sps)
        if cre1 * cre2 * ani1 * ani2 != 0
            break 
        end
        onehot = Int128(1) << (idx_m-1)
        if ani1 == 0 && (bitint_ket & onehot) == onehot && (bitint_bra & onehot) == 0
            ani1 = idx_m
        end
        if ani1 != 0 && (bitint_ket & onehot) == onehot && (bitint_bra & onehot) == 0 && idx_m != ani1
            ani2 = idx_m
        end
        if cre1 == 0 && (bitint_bra & onehot) == onehot && (bitint_ket & onehot) == 0
            cre1 = idx_m
        end
        if cre1 != 0 && (bitint_bra & onehot) == onehot && (bitint_ket & onehot) == 0 && idx_m != cre1
            cre2 = idx_m
        end
    end

    phase_bra = phase_ket = 1
    phase_from_ani1 = phase_from_ani2 = phase_from_cre1 = phase_from_cre2 = 1
    for idx_m in 1:length(m_sps)
        onehot = Int128(1) << (idx_m-1)
        if (bitint_ket & onehot) == onehot
            if idx_m < ani1; phase_ket *= -1; phase_from_ani1 *= -1 ; end
            if idx_m < ani2 && idx_m != ani1; phase_ket *= -1; phase_from_ani2 *= -1; end
        end
        if (bitint_bra & onehot) == onehot
            if idx_m < cre1; phase_bra *= -1; phase_from_cre1 *= -1;end
            if idx_m < cre2 && idx_m != cre1; phase_bra *= -1; phase_from_cre2 *= -1;  end
        end
    end
    phase = phase_bra * phase_ket 

    idx_jj_p = dict_m2j[cre1 + offset_idx]
    idx_jj_q = dict_m2j[cre2 + offset_idx]
    idx_jj_r = dict_m2j[ani1 + offset_idx]
    idx_jj_s = dict_m2j[ani2 + offset_idx]

    j1 = m_sps[cre1].j; m1 = m_sps[cre1].jz
    j2 = m_sps[cre2].j; m2 = m_sps[cre2].jz
    j3 = m_sps[ani1].j; m3 = m_sps[ani1].jz
    j4 = m_sps[ani2].j; m4 = m_sps[ani2].jz

    txt = ""; txt_count = 0 
    if verbose >= 1
        txt  = "op ~ $(cre1)_+ $(cre2)_+ | $(ani1)_- $(ani2)_- \n"
        txt *= "jm_1 ($j1, $m1) jm_2 ($j2, $m2) jm_3 ($j3, $m3) jm_4 ($j4, $m4) \n"
        txt *= "phase bra $phase_bra ket $phase_ket => $phase\n" 
        txt *= "phase_from cre1 $phase_from_cre1 cre2 $phase_from_cre2\n" 
        txt *= "phase_from ani1 $phase_from_ani1 ani2 $phase_from_ani2" 
    end
    H_ij = 0.0
    key_bra = hash_2ints(idx_jj_p, idx_jj_q)
    key_ket = hash_2ints(idx_jj_r, idx_jj_s)
    if haskey(hamil, (key_bra, key_ket))
        relevant = hamil[(key_bra, key_ket)]
        for (J, Vjj) in relevant
            J = Int(J)
            Mbra = div(m1+m2, 2); if abs(Mbra) > J; continue; end
            Mket = div(m3+m4, 2); if abs(Mket) > J; continue; end        
            NJ_pq = sqrt(deltaf(idx_jj_p, idx_jj_q) * (-1)^J + 1.0) / (1 + deltaf(idx_jj_p, idx_jj_q))
            NJ_rs = sqrt(deltaf(idx_jj_r, idx_jj_s) * (-1)^J + 1.0) / (1 + deltaf(idx_jj_r, idx_jj_s))
            cg1 = dict_CGs[get_nkey6_shift(j1, m1, j2, m2, J, Mbra; int_shift=int_shift)]
            cg2 = dict_CGs[get_nkey6_shift(j3, m3, j4, m4, J, Mket; int_shift=int_shift)]
            CG = cg1 * cg2
            vtmp = Vjj * CG / (NJ_pq * NJ_rs) *phase
            H_ij += vtmp
            if verbose >= 1
                txt_count += 1
                if txt_count == 1
                    println(txt)
                end
                println("Vjj(04/40): <$idx_jj_p $idx_jj_q| $idx_jj_r $idx_jj_s> J=$J v=$Vjj v*fac=$(CG*Vjj/(NJ_pq*NJ_rs)) cg1 $cg1 cg2 $cg2 NJ_pq=$NJ_pq NJ_rs=$NJ_rs")
            end
        end
    end
    return H_ij
    
end

function eval_Hij_diag(Hamil_snt, p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket,
                       p_msps, n_msps, dict_m2j, int_shift, dict_CGs; verbose::Bool=false)
    H_ij = H_1b = H_pp = H_nn = H_pn = 0.0
    # proton/neutron 1-body term
    for idx_m in 1:length(p_msps)
        idx_j = dict_m2j[idx_m]
        onehot = Int128(1) << (idx_m-1)
        if (p_bitint_bra & onehot) == onehot == (p_bitint_ket & onehot) 
            H_1b += Hamil_snt.h_1b[(idx_j, idx_j)]
        end
    end
    offset = length(p_msps)
    for idx_m in 1:length(n_msps)
        idx_j = dict_m2j[idx_m+offset]
        onehot = Int128(1) << (idx_m-1)
        if (n_bitint_bra & onehot) == onehot == (n_bitint_ket & onehot) 
            H_1b += Hamil_snt.h_1b[(idx_j, idx_j)]
        end
    end
    H_ij += H_1b

    # Vpp&Vnn terms
    for ppnn = 1:2
        msps = ifelse(ppnn == 1, p_msps, n_msps)
        bitint_bra = ifelse(ppnn == 1, p_bitint_bra, n_bitint_bra)
        bitint_ket = ifelse(ppnn == 1, p_bitint_ket, n_bitint_ket)
        ofst = ifelse(ppnn == 1, 0, length(p_msps))
        for idx_m_1 in 1:length(msps)
            for idx_m_2 in idx_m_1+1:length(msps)
                idx_j_1 = dict_m2j[idx_m_1 + ofst]
                idx_j_2 = dict_m2j[idx_m_2 + ofst]
                onehot_1 = Int128(1) << (idx_m_1-1)
                onehot_2 = Int128(1) << (idx_m_2-1)
                if  (bitint_bra & onehot_1) != onehot_1 ||
                    (bitint_bra & onehot_2) != onehot_2 ||
                    (bitint_ket & onehot_1) != onehot_1 ||
                    (bitint_ket & onehot_2) != onehot_2
                    continue
                end
                target = ifelse(ppnn == 1, Hamil_snt.V2b_pp, Hamil_snt.V2b_nn)
                # println("idx_i_m $idx_m_1 idx_i_j $idx_j_1 ; idx_j_m $idx_m_2  idx_j_j $idx_j_2")
                # println("target keys: $(keys(target))")
                key1 = hash_2ints(idx_j_1, idx_j_2)
                key2 = hash_2ints(idx_j_1, idx_j_2)
                if !haskey(target, (key1, key2))
                    continue
                end
                Vs_relevant = target[(key1, key2)]
                j1 = msps[idx_m_1].j; m1 = msps[idx_m_1].jz
                j2 = msps[idx_m_2].j; m2 = msps[idx_m_2].jz
                M = div(m1+m2, 2)
                for (J, V_jj) in Vs_relevant
                    J = Int(J)
                    Nfactor = sqrt(deltaf(idx_j_1, idx_j_2) * (-1)^J + 1.0) / (1 + deltaf(idx_j_1, idx_j_2))
                    if abs(div(m1+m2,2)) > J
                        continue
                    end                    
                    nkey = get_nkey6_shift(j1, m1, j2, m2, J, M; int_shift=int_shift)
                    CG = dict_CGs[nkey]
                    vtmp = V_jj * (CG / Nfactor)^2  
                    if verbose
                        println("Vjj $V_jj J $J vtmp $vtmp CG $CG ")
                    end
                    H_ij += vtmp
                    if ppnn == 1
                        H_pp += vtmp
                    else
                        H_nn += vtmp
                    end
                end
            end
        end
    end
    if verbose
        println("H_ij = $H_ij, H_1b = $H_1b, H_pp = $H_pp, H_nn = $H_nn")
    end
    return H_ij 
end

function eval_Hij_pn(Hamil_snt, p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket, 
                     p_ham_dist, n_ham_dist, p_msps, n_msps, dict_m2j, int_shift, dict_CGs, verbose, to)
    hamil = Hamil_snt.V2b_pn
    H_ij = 0.0

    cre_p = ani_p = cre_n = ani_n = 0
    phase_bra_p = phase_ket_p = phase_bra_n = phase_ket_n = 1
    #@timeit to "pn_idxloop" 
    for pn_idx in 1:2
        m_sps = ifelse(pn_idx == 1, p_msps, n_msps)
        bitint_bra = ifelse(pn_idx == 1, p_bitint_bra, n_bitint_bra)
        bitint_ket = ifelse(pn_idx == 1, p_bitint_ket, n_bitint_ket)
        for idx_m in 1:length(m_sps)
            onehot = Int128(1) << (idx_m-1)
            if (bitint_ket & onehot) == onehot && (bitint_bra & onehot) == 0
                if pn_idx == 1
                    ani_p = idx_m 
                else
                    ani_n = idx_m 
                end
            end
            if (bitint_bra & onehot) == onehot && (bitint_ket & onehot) == 0
                if pn_idx == 1
                    cre_p = idx_m 
                else
                    cre_n = idx_m 
                end
            end
        end
        for idx_m in 1:length(m_sps) 
            onehot = Int128(1) << (idx_m-1)
            if (bitint_ket & onehot) == onehot 
                if pn_idx == 1
                    if idx_m < ani_p 
                        phase_ket_p *= -1
                    end                    
                else
                    if idx_m < ani_n
                        phase_ket_n *= -1
                    end
                end
            end
            if (bitint_bra & onehot) == onehot
                if pn_idx == 1
                    if idx_m < cre_p                        
                        phase_bra_p *= -1
                    end
                else
                    if idx_m < cre_n
                        phase_bra_n *= -1
                    end
                end
            end
        end
    end
    ps = [ (cre_p, ani_p)]
    ns = [ (cre_n, ani_n)]
    if cre_p == ani_p == 0
        ps = [ (i, i) for i in 1:length(p_msps) if ( p_bitint_bra & (Int128(1) << (i-1)) ) == Int128(1) << (i-1) ]
    end
    if cre_n == ani_n == 0
        ns = [ (i, i) for i in 1:length(n_msps) if ( n_bitint_bra & (Int128(1) << (i-1)) ) == Int128(1) << (i-1) ]
    end
    if verbose >= 1
        println("cre_p $cre_p ani_p $ani_p cre_n $cre_n ani_n $ani_n ",
                "phase_bra_p $phase_bra_p phase_ket_p $phase_ket_p phase_bra_n $phase_bra_n phase_ket_n $phase_ket_n")
    end
    phase_fac = phase_bra_p * phase_ket_n * phase_bra_n * phase_ket_p
    for (cre_p, ani_p) in ps
        idx_jj_ani_p = dict_m2j[ani_p]
        idx_jj_cre_p = dict_m2j[cre_p]
        j_p_bra = p_msps[cre_p].j; m_p_bra = p_msps[cre_p].jz
        j_p_ket = p_msps[ani_p].j; m_p_ket = p_msps[ani_p].jz
        key_p = hash_2ints(idx_jj_cre_p, idx_jj_ani_p)
        for (cre_n, ani_n) in ns
            idx_jj_ani_n = dict_m2j[ani_n + length(p_msps)]
            idx_jj_cre_n = dict_m2j[cre_n + length(p_msps)]
            j_n_bra = n_msps[cre_n].j; m_n_bra = n_msps[cre_n].jz
            j_n_ket = n_msps[ani_n].j; m_n_ket = n_msps[ani_n].jz
            key_n = hash_2ints(idx_jj_cre_n, idx_jj_ani_n)
            Mbra = div(m_p_bra + m_n_bra, 2)
            Mket = div(m_p_ket + m_n_ket, 2)
            if verbose >= 1
                println("op~ $(cre_p)_+ $(cre_n+length(p_msps))_+ $(ani_p)_- $(ani_n+length(p_msps))_- ",
                         "in jj : $idx_jj_cre_p $idx_jj_cre_n $idx_jj_ani_p $idx_jj_ani_n")
            end
            if haskey(hamil, (key_p, key_n))
                relevant = hamil[(key_p, key_n)]
                for (J, Vpn) in relevant
                    J = Int(J) 
                    if abs(Mbra) > J || abs(Mket) > J
                        continue
                    end
                    CG  = dict_CGs[get_nkey6_shift(j_p_bra, m_p_bra, j_n_bra, m_n_bra, J, Mbra; int_shift=int_shift)]
                    CG *= dict_CGs[get_nkey6_shift(j_p_ket, m_p_ket, j_n_ket, m_n_ket, J, Mket; int_shift=int_shift)]
                    vtmp = Vpn * CG  * phase_fac
                    H_ij += vtmp
                    if verbose >= 1
                        println("Vpn: J=$J v=$(@sprintf("%10.6f",Vpn)) vtmp=$(@sprintf("%10.6f", vtmp)) CG $(@sprintf("%10.6f",CG)) ")
                    end
                end
            end
        end
    end
    return H_ij
end

"""
    operate_H_on_vec!(w, Hamil_mat::Array{Float64, 2}, v)

Function to compute the matrix-vector product of the Hamiltonian matrix and a vector.
"""
function operate_H_on_vec!(w, Hamil_mat::Array{Float64, 2}, v)
    mul!(w, Hamil_mat, v)
    return nothing
end

"""
    operate_H_on_vec!(w, Hamil_mat::Dict{UInt64, Float64}, v)

Sparse version of the function to compute the matrix-vector product of the Hamiltonian matrix in the form of `Dict{UInt64, Float64}` and a vector.
"""
function operate_H_on_vec!(w, Hamil_mat::Dict{UInt64, Float64}, v, partials::Matrix{Float64}, keysvec, to)
    Threads.@threads for k in eachindex(keysvec)
        tid = Threads.threadid()
        partial = @view partials[:, tid]
        nkey = keysvec[k]
        val = Hamil_mat[nkey]
        i, j = unhash_2ints(nkey)
        partial[i] += val * v[j]
        if i != j
            partial[j] += val * v[i]
        end
    end
    # Sum the thread-local accumulators into the output vector.
    w .= 0.0
    @inbounds for tid in 1:size(partials, 2)
        w .+= @view partials[:, tid]
    end
    partials .= 0.0
    return nothing
end

"""
    lanczos(Hamil_mat, dim, get_evecs, to; itnum=300, tol=1e-9, debug_mode=0)

Function to compute the lowest eigenvalue of the Hamiltonian using the Lanczos method.

Constructing a Krylov subspace ``\\mathcal{K}_m(H,v) = \\mathrm{span}\\{v, Hv, H^2v, \\cdots, H^{m-1}v\\}``,
and the tridiagonal matrix ``T_m = V_m^T H V_m`` where ``V_m = [v_1, v_2, \\cdots, v_m]`` and ``v_{m+1} = H v_m - \\alpha_m v_m - \\beta_m v_{m-1}``, the Lanczos method iteratively constructs the matrix ``T_m`` and diagonalizes it to obtain the smallest eigenvalue of ``H``.

# Arguments
- `Hamil_mat`: Hamiltonian matrix, either a `Matrix{Float64}` or a `Dict{UInt64, Float64}`.
- `dim`: Dimension of the basis, i.e. number of configurations
- `get_evecs`: If `true`, save the exact wave function to a HDF5 file
- `to`: TimerOutput object

# Optional arguments
- `itnum`: Number of Lanczos iterations
- `tol`: Tolerance for convergence
- `debug_mode`: the level of debug information
"""
function lanczos(Hamil_mat, dim, n_eigen, get_evecs, to; itnum=300, tol=1e-8, debug_mode=0)
    if n_eigen >= dim
        n_eigen = dim
    end

    keysvec = collect(keys(Hamil_mat))
    partials = zeros(Float64, dim, Threads.maxthreadid())
    println("Starting Lanczos iteration...")
    Random.seed!(123456)
    v1 = rand(dim)
    v1 ./= norm(v1)
    alpha = beta = 0.0
    T = zeros(Float64, itnum, itnum)
    w = zeros(Float64, dim)
    Es_monitor = zeros(Float64, n_eigen, 2)

    vks = [zeros(Float64, dim) for _ in 1:itnum]
    vks[2] = v1
    conv_flag = 0
    it_finished = min(dim, itnum)
    for i = 1:min(dim, itnum)
        conv_flag = 0
        v = vks[i+1]
        @timeit to "operate" operate_H_on_vec!(w, Hamil_mat, v, partials, keysvec, to)
        alpha = dot(w, v)
        T[i, i] = alpha
        if i >= n_eigen
            Es_monitor[:, 1] .= eigvals(T[1:i, 1:i])[1:n_eigen]
            if (debug_mode > 0 && i % 10 == 0) || (debug_mode > 1)
                print_vec("iter = $(@sprintf("%6i", i))", Es_monitor[:,1])
            end
            for n in 1:n_eigen
                if abs(Es_monitor[n, 1] - Es_monitor[n, 2]) < tol
                    conv_flag += 1
                end
            end
        end
        it_finished = i
        if conv_flag == n_eigen
            break
        end
        w .-= alpha .* v
        w .-= beta .* vks[i]
        @timeit to "reOrth" reOrthogonalize!(w, vks, i)
        beta = norm(w)
        vks[i+2] = w ./ beta
        T[i, i+1] = beta
        T[i+1, i] = beta
        if i >= n_eigen
            Es_monitor[:, 2] .= Es_monitor[:, 1]
        end
    end
    evecs = zeros(Float64, 1, 1)
    if get_evecs
        evecs = get_ritz_vector(vks, T, n_eigen, it_finished)
    end
    return Es_monitor[:, 1], evecs
end

function get_ritz_vector(vks, Tmat, n_eigen, it_finished)
    dim = length(vks[1])
    vals, vecs = eigen(@views Tmat[1:it_finished, 1:it_finished])
    Rvecs = zeros(Float64, dim, n_eigen)
    for nth in 1:n_eigen
        Rvec = @view Rvecs[:, nth]
        for k in 1:size(vecs, 1)
            coeff = vecs[k, nth]
            if k+1 > length(vks)
                break
            end
            axpy!(coeff, vks[k+1], Rvec)
        end
        Rvec .*= 1.0/sqrt(dot(Rvec,Rvec))
    end    
    return Rvecs
end

"""
    reOrthogonalize!(w, vks, i)

Re-orthogonalize the vector w with respect to the previous vectors:

```math
w := w - \\sum_{j=1}^{i} \\langle w, v_j \\rangle v_j
```
"""
function reOrthogonalize!(w, vks, i)
    for j in 1:i
        w .-= dot(w, vks[j]) .* vks[j]
    end
    return w
end


function prepare_CGcoeffs(p_msps, n_msps)
    dict_CGs = Dict{UInt, Float64}( )
    jmax = maximum( p_msps[n].j for n in 1:length(p_msps) )
    jmax = max(jmax, maximum( n_msps[n].j for n in 1:length(n_msps) ))
    for j1 = 1:2:jmax
        for j2 = 1:2:jmax
            Jmin = abs(j1 - j2)
            Jmax = j1 + j2
            for J in Jmin:2:Jmax
                J = div(J, 2)
                for m1 in -j1:2:j1
                    for m2 in -j2:2:j2
                        M = div(m1 + m2, 2)
                        if abs(M) > J; continue; end
                        CG = clebschgordan(Float64, j1//2, m1//2, j2//2, m2//2, J//1, M//1)
                        key = get_nkey6_shift(j1, m1, j2, m2, J, M; int_shift=jmax+1)
                        dict_CGs[key] = CG
                    end
                end
            end
        end
    end
    println("CG coefficients prepared (jmax=$(jmax)), total $(length(dict_CGs)) entries.")
    return jmax+1, dict_CGs
end


function construct_Hmat(Hamil_snt, mdim, all_bitint_prod, p_msps, n_msps, vZ, vN, 
                        dict_m2j, int_shift, dict_CGs, verbose, to;
                        Hrank=2)

    Hmat = Dict{UInt, Float64}( )
    partial_dict = [ Dict{UInt, Float64}( ) for _ in 1:Threads.maxthreadid() ]

    # # 3NF stuff
    is_3NF = length(Hamil_snt.V3b) > 0
    #Hamil_3NF = Hamil_snt.V3b[1]
    # println("checking 3NF sps...")
    # sps_3b = Hamil_3NF.sps_3b.sps
    # for (idx, sp) in sps_3b
    #     println("3NF sps idx $idx : n $(sp.n) l $(sp.l) j $(sp.j) tz $(sp.tz)")
    # end
    # println("checking 3NF terms...")
    # cnt = 0
    # for (key, v3m) in Hamil_3NF.v3monopole        
    #     a, b, c, d, e, f = unhash_key6j(key)
    #     println("<$a $b $c |V3b| $d $e $f > = $v3m")
    #     cnt += 1
    #     if cnt == 5
    #         break
    #     end
    # end

    @timeit to "thread loop" Threads.@threads for iter in 1:mdim^2
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
        prod_ham_dist = p_ham_dist * n_ham_dist
        if sum_ham_dist > Hrank * 2 # bra and ket cannot be connected by the Hamiltonian
            continue 
        end

        tid = Threads.threadid()
        # # Calculate the Hamiltonian matrix element
        Hij = 0.0

        # sum_ham_dist = 0 case, i.e. bra and ket are the same
        #@timeit to "diag" 
        if sum_ham_dist == 0 
            Hij += eval_Hij_diag(Hamil_snt, p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket,
                                 p_msps, n_msps, dict_m2j, int_shift, dict_CGs)
        end

        #@timeit to "02_20" 
        if sum_ham_dist == 2
            Hij += eval_Hij_nondiag_02_20(Hamil_snt, p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket,
                                          p_ham_dist, n_ham_dist, p_msps, n_msps, dict_m2j, Hrank, int_shift, dict_CGs, verbose)
        end 

        #@timeit to "04_40" 
        if sum_ham_dist == 4 && prod_ham_dist == 0
            Hij += eval_Hij_nondiag_04_40(Hamil_snt, p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket,
                                          p_ham_dist, n_ham_dist, p_msps, n_msps, dict_m2j, Hrank, int_shift, dict_CGs,verbose)
        end

        #@timeit to "pn" 
        if vZ * vN != 0 && ( sum_ham_dist <= 2 || (sum_ham_dist == prod_ham_dist == 4) )
            Hij += eval_Hij_pn(Hamil_snt, p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket, 
                                p_ham_dist, n_ham_dist, p_msps, n_msps, dict_m2j, int_shift, dict_CGs, verbose, to) 
        end
        
        # Dev not: I am going to add 3NF here
        # First, need to evaluate T=3/2 terms (ppp/nnn),
        # then T=1/2 terms (ppn/nnp) will be considered... long way to go.

        if is_3NF
            #@timeit to "3NF T=3/2"
            Hij += eval_Hij_3NF_T3(Hamil_snt, p_bitint_bra, n_bitint_bra, p_bitint_ket, n_bitint_ket,
                                   p_ham_dist, n_ham_dist, p_msps, n_msps, dict_m2j, Hrank, int_shift, dict_CGs, verbose)
        end


        #@timeit to "assign" 
        partial_dict[tid][hash_2ints(idx_bra, idx_ket)] = Hij
    end    

    # Combine the partial results from all threads
    @timeit to "collect" for tid in 1:Threads.maxthreadid()
        for (nkey, val) in partial_dict[tid]
            Hmat[nkey] = val
        end
    end
    if verbose >=1
        println("Hamiltonian matrix")
        for (nkey, val) in Hmat
            idx_bra, idx_ket = unhash_2ints(nkey)
            bra_bitint, ket_bitint = all_bitint_prod[idx_bra], all_bitint_prod[idx_ket]
            p_bra, n_bra = bra_bitint
            p_ket, n_ket = ket_bitint
            bra_bitstr = int2bitstr(p_bra, length(p_msps)) * " ⊗ " * int2bitstr(n_bra, length(n_msps))
            ket_bitstr = int2bitstr(p_ket, length(p_msps)) * " ⊗ " * int2bitstr(n_ket, length(n_msps))
            println("H[$(idx_bra),$(idx_ket)] = $(Hmat[nkey]) for bra $(bra_bitstr) ket $(ket_bitstr)")

        end
    end
    return Hmat
end

function prepare_all_configs_in_modelspace(parity, Mtot, vZ, vN, p_msps, n_msps, Nref_proton, Nref_neutron, verbose, to)
    @timeit to "configs." begin
        struct_proton_config  = count_only_num_particle(vZ, p_msps, Nref_proton; store_config=true)
        struct_neutron_config = count_only_num_particle(vN, n_msps, Nref_neutron; store_config=true)
        if verbose >= 2
            for (i, MPNblock) in enumerate(struct_proton_config)
                println("Proton $(i)-th block: Mtot = $(MPNblock.Mtot), parity = $(MPNblock.parity) Nexc = $(MPNblock.Nexc) Dim. = $(MPNblock.count[1])")
                println("configs: $(MPNblock.configs)")
                @assert length(MPNblock.configs) == MPNblock.count[1] "Count mismatch for proton block $(i)"
            end
            for (i, MNNblock) in enumerate(struct_neutron_config)
                println("Neutron $(i)-th block: Mtot = $(MNNblock.Mtot), parity = $(MNNblock.parity) Nexc = $(MNNblock.Nexc) Dim. = $(MNNblock.count[1])")
                println("configs: $(MNNblock.configs)")
                @assert length(MNNblock.configs) == MNNblock.count[1] "Count mismatch for neutron block $(i)"
            end
        end
    end
    num_block_p = length(struct_proton_config)
    num_block_n = length(struct_neutron_config)

    all_bitint_prod = Tuple{Int, Int}[ ]
    mdim = 0
    @timeit to "Counting M-dim." for bidx_p in 1:num_block_p
        p_block = struct_proton_config[bidx_p]
        parity_proton = p_block.parity
        Mp = p_block.Mtot
        for bidx_n in 1:num_block_n
            n_block = struct_neutron_config[bidx_n]
            parity_neutron = n_block.parity
            if parity_proton * parity_neutron != parity
                continue
            end
            Mn = n_block.Mtot
            if Mp + Mn != Mtot
                continue
            end
            tdim = p_block.count[1] * n_block.count[1]
            mdim += tdim

            for pbit in p_block.configs
                for nbit in n_block.configs
                    push!(all_bitint_prod, (pbit, nbit))
                end
            end
        end      
    end
    mdim = Int(mdim)
    print("Full M-scheme dimension: $(mdim)    ")
    
    if verbose >= 1 
        println("Checking the bitstrings and their properties:")
        for (i, bitint_prod) in enumerate(all_bitint_prod)
            p_bitint, n_bitint = bitint_prod
            mconfig_p = read_bitint(p_bitint, p_msps; Nexc_ref=Nref_proton)
            mconfig_n = read_bitint(n_bitint, n_msps; Nexc_ref=Nref_neutron)
            Nocc_p = mconfig_p.Nocc; Nocc_n = mconfig_n.Nocc
            Nocc = Nocc_p + Nocc_n
            Nexc_p = mconfig_p.Nexc; Nexc_n = mconfig_n.Nexc
            Nexc = Nexc_p + Nexc_n
            Mp = mconfig_p.Mtot; Mn = mconfig_n.Mtot
            tMtot = Mp + Mn
            tparity = mconfig_p.parity * mconfig_n.parity

            @assert Nocc == Nocc_p + Nocc_n "Nocc mismatch: $Nocc != $Nocc_p + $Nocc_n"
            @assert Mtot == tMtot "Mtot mismatch: $Mtot != $tMtot"
            @assert parity == tparity "Parity mismatch: $parity != $(tparity)"
            if verbose >= 3
                println(@sprintf("%10i", i), ": ", bitstr, " => Nocc = $(Nocc) (p: $(Nocc_p), n: $(Nocc_n)), Mtot = $(Mtot) (Mp $(Mp), Mn $(Mn)), parity = $(parity), Nexc = $(Nexc)")
            end
        end
    end
    return all_bitint_prod
end

function summarize_bitstrings(sampled_bits, p_msps, n_msps, vZ, vN, maxnum_subspace_basis::Int;
                              Qiskit_ordered=true, Nref_proton=0, Nref_neutron=0, 
                              postselection::Vector{Int}=[0, 0, 0],
                              verbose::Bool = false,
                              preferred_configs::Vector{Tuple{Int, Int}} = [ ])
    Nocc, parity, Mtot = postselection
    sampled_bits = sampled_bits[1:min(maxnum_subspace_basis, length(sampled_bits))]
    all_bitint_prod = Tuple{Int128, Int128}[ ]
    if length(preferred_configs) > 0
        all_bitint_prod = [ (Int128(p_bit), Int128(n_bit)) for (p_bit, n_bit) in preferred_configs ]
    end
    if typeof(sampled_bits) == Vector{Int} || typeof(sampled_bits) == Vector{Int128}
        @assert Qiskit_ordered "sample_bits::Vector{Int} or Vector{Int128} is not supported in Qiskit_ordered=false mode."
        for bit in sampled_bits
            pbitint = Int128(bit) >> length(n_msps)
            nbitint = Int128(bit) & ((Int128(1) << length(n_msps)) - 1)
            if (pbitint, nbitint) in all_bitint_prod
                continue
            end
            push!(all_bitint_prod, (pbitint, nbitint))
            msg = "The number of particles in the bitstring does not match the expected vZ and vN. "
            msg *= "pbit $(int2bitstr(pbitint, length(p_msps))) nbit $(int2bitstr(nbitint, length(n_msps)))"
            @assert count_ones(pbitint) == vZ && count_ones(nbitint) == vN msg
        end
    elseif typeof(sampled_bits) == Vector{String}
        @assert length(sampled_bits[1]) == length(p_msps) + length(n_msps) "Each bitstring must have length $(length(p_msps) + length(n_msps)), but got $(length(sampled_bits[1]))"
        targets = nothing
        if Qiskit_ordered # neutron(<=) ⊗ proton(<=)
            targets = [(parse(Int128, bit[length(n_msps)+1:end], base=2),
                        parse(Int128, bit[1:length(n_msps)], base=2)) for bit in sampled_bits]
        else # neutron(=>) ⊗ proton(=>) this may not be used 
            targets = [(parse(Int128, bit[end:-1:length(n_msps)+1], base=2),
                        parse(Int128, bit[length(n_msps):-1:1], base=2)) for bit in sampled_bits]
        end
        for bitint_prod in targets
            pbitint, nbitint = bitint_prod
            if (pbitint, nbitint) in all_bitint_prod
                continue
            end
            push!(all_bitint_prod, (pbitint, nbitint))
            msg = "The number of particles in the bitstring does not match the expected vZ and vN. "
            msg *= "pbit $(int2bitstr(pbitint, length(p_msps))) nbit $(int2bitstr(nbitint, length(n_msps)))"
            @assert count_ones(pbitint) == vZ && count_ones(nbitint) == vN msg
        end
    end

    if postselection != [0, 0, 0]
        truncated = Tuple{Int128, Int128}[ ]
        for (i, bitint_prod) in enumerate(all_bitint_prod)
            p_bitint, n_bitint = bitint_prod
            mconfig_p = read_bitint(p_bitint, p_msps; Nexc_ref=Nref_proton)
            mconfig_n = read_bitint(n_bitint, n_msps; Nexc_ref=Nref_neutron)
            Nocc_p = mconfig_p.Nocc; Nocc_n = mconfig_n.Nocc
            tNocc = Nocc_p + Nocc_n
            Nexc_p = mconfig_p.Nexc; Nexc_n = mconfig_n.Nexc
            Nexc = Nexc_p + Nexc_n
            Mp = mconfig_p.Mtot; Mn = mconfig_n.Mtot
            tMtot = Mp + Mn
            tparity = mconfig_p.parity * mconfig_n.parity
            if verbose
                println("$(int2bitstr(p_bitint, length(p_msps))) ⊗ $(int2bitstr(n_bitint, length(n_msps))) => ",
                        "Nocc = $(tNocc) (p: $(Nocc_p), n: $(Nocc_n)), Mtot = $(tMtot) (Mp $(Mp), Mn $(Mn)), parity = $(tparity) ")
            end
            if tNocc != Nocc; continue; end
            if tMtot != Mtot; continue; end
            if tparity != parity; continue; end
            push!(truncated, bitint_prod)
        end
        all_bitint_prod = truncated
        truncated = nothing
    end

    if verbose
        println("checking sampled bitstrings: M=$(Mtot), parity=$(parity), Nocc=$(Nocc)")
        if length(all_bitint_prod) == length(sampled_bits)
            for (i, bit_input) in enumerate(sampled_bits)
                println("bitinput $bit_input")
                println("$bit_input => $(int2bitstr(all_bitint_prod[i][1], length(p_msps))) ⊗ $(int2bitstr(all_bitint_prod[i][2], length(n_msps)))")
            end
        else
            for (pbit, nbit) in all_bitint_prod
                println("$(int2bitstr(pbit, length(p_msps))) ⊗ $(int2bitstr(nbit, length(n_msps)))")
            end
        end
    end
    return all_bitint_prod
end

struct Result_QSCI
    evals::Vector{Float64}
    evecs::Matrix{Float64}
    Jvals::Vector{Float64}
    evars::Vector{Float64}
    mdim::Int
end

function evaluate_energy_variance(Hamil_mat, evecs, n_eigen, to)
    keysvec = collect(keys(Hamil_mat))
    evars = zeros(Float64, n_eigen)
    dim = size(evecs, 1)
    ws = zeros(Float64, dim, n_eigen)
    us = zeros(Float64, dim, n_eigen)
    partials = zeros(Float64, dim, Threads.maxthreadid())
    for i in 1:n_eigen
        w = @view ws[:, i]
        u = @view us[:, i]
        v = @view evecs[:, i]
        operate_H_on_vec!(u, Hamil_mat, v, partials, keysvec, to)       
        operate_H_on_vec!(w, Hamil_mat, u, partials, keysvec, to)
        evars[i] = dot(w, v)
    end
    return evars
end

function get_occs_jj(SPEs, sps_jj, Nocc, NpNh=0)
    @assert 0 <= NpNh <= 2 "NpNh must be 0, 1, or 2 in the current implementation."
    len = length(SPEs)
    argmins = sortperm(SPEs)
    occs_jj = zeros(Int, len)
    org_idx_Fermi_level = idx_Fermi_level = count = 0
    for idx in 1:len
        idx_jj = argmins[idx]
        j2 = sps_jj[idx_jj].j        
        if count + j2 + 1 >= Nocc
            idx_Fermi_level = idx 
            org_idx_Fermi_level = argmins[idx_Fermi_level]
        end        
        #println("idx: $idx j2: $j2 count: $count Nocc: $Nocc occs_jj: $occs_jj")
        for jz = -j2:2:j2
            count += 1
            occs_jj[idx_jj] += 1
            if count == Nocc
                break
            end
        end
        if count == Nocc
            break
        end
    end
    occs_summary = Vector{Int}[ ]
    push!(occs_summary, occs_jj)
    if NpNh > 0 
        occs_summary = Vector{Int}[ ]
        for idx_particle in argmins[idx_Fermi_level+1:end]
            occs = zeros(Int, len) .+ occs_jj
            oidx_particle = argmins[idx_particle]
            #println("argmins $argmins idx_Fermi_level (@$(idx_Fermi_level)): $org_idx_Fermi_level idx_particle (@$(idx_particle)): $oidx_particle") 
            if sps_jj[idx_particle].j + 1 >= NpNh && sps_jj[org_idx_Fermi_level].j + 1 >= NpNh
                occs[org_idx_Fermi_level] -= NpNh
                occs[idx_particle] += NpNh
            else
                @error "This case NpNh=$NpNh is not supported."
            end
            push!(occs_summary, occs)
        end
    end
    return occs_summary
end

function get_lowest_filling_configs(SPEs, sps_jj, Nocc, dict_j2m, include_2p2h::Bool)
    len = length(SPEs)

    occs_pool = get_occs_jj(SPEs, sps_jj, Nocc, 0)
    if include_2p2h
        occs_pool = vcat(occs_pool, get_occs_jj(SPEs, sps_jj, Nocc, 1))
        occs_pool = vcat(occs_pool, get_occs_jj(SPEs, sps_jj, Nocc, 2))
    end
    occ_full = [ sps_jj[i].j + 1 for i in 1:len ]
    len_m = sum( sps_jj[i].j + 1 for i in 1:len )
    bitint_possible = Int128[ ]
    for occs_jj in occs_pool
        bitint_base = Int128(0)
        for idx_jj in 1:len
            if occs_jj[idx_jj] == 0
                continue
            end
            # jj is fully occupied case
            if occs_jj[idx_jj] == sps_jj[idx_jj].j +1 
                for idx_mm in dict_j2m[idx_jj]
                    bitint_base = bitint_base | (Int128(1) << (idx_mm - 1))
                end
            end
        end 
        #println("occs_jj: $occs_jj, bitint_base: $bitint_base $(int2bitstr(bitint_base, len_m))")
        count_partial_jj = 0
        subbit_pool = [ Int128[ ] for _ in 1:len ]
        for idx_jj in 1:len
            if occs_jj[idx_jj] == 0
                continue 
            end
            Nocc_jj = occs_jj[idx_jj]
            # partially occupied case
            if occs_jj[idx_jj] < sps_jj[idx_jj].j + 1
                count_partial_jj += 1
                x = bitint_base
                idx_min_relevant = minimum(dict_j2m[idx_jj])
                bit_rel_min = Int128(0)
                nshift = sps_jj[idx_jj].j + 1 - Nocc_jj
                #println("nshift: $nshift  idx_min_relevant: $idx_min_relevant ")
                for idx_mm = idx_min_relevant: idx_min_relevant + Nocc_jj - 1
                    bit_rel_min |= (Int128(1) << (idx_mm - 1))
                end
                bit_rel_max = bit_rel_min << nshift
                #println("bit_rel_min: $bit_rel_min $(int2bitstr(bit_rel_min, len_m))")
                #println("bit_rel_max: $bit_rel_max $(int2bitstr(bit_rel_max, len_m))")

                nshift_sub = idx_min_relevant - 1
                x = bit_rel_min >> nshift_sub # Integer within target jj orbit
                push!(subbit_pool[idx_jj], x << nshift_sub | bitint_base)
                while x < (bit_rel_max >> nshift_sub)
                    u = x & (-x)
                    v = x + u   
                    x = v + (((v ⊻ x) ÷ u) >> 2) 
                    xtarget = x << nshift_sub | bitint_base
                    #println("xtarget: $xtarget $(int2bitstr(xtarget, len_m))")
                    push!(subbit_pool[idx_jj], xtarget)
                end
            end
        end
        if count_partial_jj == 0
            push!(bitint_possible, bitint_base)
        end

        ## Here, we have to sum up all the possible combinations of subbit_pool
        ## For example, if subbit_pool = [ [2, 32], [4, 16, 64], [], [ ] ], then we have 2+4, 2+16, 2+64, 32+4, 32+16, 32+64.
        ## i.e. picking one from each subbit_pool
        #println("subbit_pool: $subbit_pool")
        if count_partial_jj > 0
            bitint_combinations = [Int128(0)]
            for idx_jj in 1:len
                if length(subbit_pool[idx_jj]) == 0
                    continue
                end
                new_combinations = [ ]
                for bitint in subbit_pool[idx_jj]
                    for existing in bitint_combinations
                        push!(new_combinations, existing | bitint)
                    end
                end
                bitint_combinations = new_combinations
            end
            for bitint in bitint_combinations
                push!(bitint_possible, bitint)
            end
        end
        #println("len(bitint_possible) = $(length(bitint_possible))")
        # for bitint in bitint_possible
        #     println("possible: $(bitint) $(int2bitstr(bitint, len_m))")
        # end
    end
    return bitint_possible
end

function random_sampling_of_configs(Hamil_snt::Hamiltonian_snt_fmt, dict_j2m,
                                    vZ, vN, all_bitint_prod, maxnum_subspace_basis, verbose;
                                    sampling_scheme::Symbol=:uniform, only_2p2h::Bool=false,
                                    nu_max::Int=10000)
    if sampling_scheme == :uniform
        idxs_subspace = sample(1:length(all_bitint_prod), maxnum_subspace_basis, replace=false)
        return idxs_subspace
    end
    if sampling_scheme == :seniority
        idxs_subspace = Int[ ]
        Dict_seniority = Dict{Int, Vector{Int}}( )
        p_msps = sps2msps(Hamil_snt.p_sps)
        n_msps = sps2msps(Hamil_snt.n_sps)
        for idx in 1:length(all_bitint_prod)
            bitint = all_bitint_prod[idx]
            vp, vn = count_seniority(bitint, p_msps, n_msps)
            if vp + vn > nu_max
                continue
            end
            if vp+vn in keys(Dict_seniority)
                push!(Dict_seniority[vp+vn], idx)
            else
                Dict_seniority[vp+vn] = [idx]
            end
        end
        keys_seniority = sort(collect(keys(Dict_seniority)))
        for key in keys_seniority
            v = key
            if verbose >= 1
                println("v = $v idx range = $(length(idxs_subspace)+1):$(length(idxs_subspace) + length(Dict_seniority[key]))")
            end
            idxs_subspace = vcat(idxs_subspace, Dict_seniority[key])
        end
        return idxs_subspace
    end
    if sampling_scheme == :lowest_filling || sampling_scheme == :lowest_filling_2p2h
        include_2p2h = (sampling_scheme == :lowest_filling_2p2h)
        p_msps = sps2msps(Hamil_snt.p_sps)
        n_msps = sps2msps(Hamil_snt.n_sps)
        # Firstly, we need to find the lowest filling configurations
        h_1b = Hamil_snt.h_1b
        lp = Hamil_snt.lp; ln = Hamil_snt.ln
        p_SPEs = [ h_1b[(i,i)] for i in 1:lp]
        n_SPEs = [ h_1b[(i+lp, i+lp)] for i in 1:ln]
        candidates_p = get_lowest_filling_configs(p_SPEs, Hamil_snt.p_sps, vZ, dict_j2m, include_2p2h)
        candidates_n = get_lowest_filling_configs(n_SPEs, Hamil_snt.n_sps, vN, dict_j2m, include_2p2h)
        if verbose >= 1
            println("candidates_p $candidates_p candidates_n $candidates_n")
        end

        idxs_subspace = Int[ ]
        for idx in 1:length(all_bitint_prod)
            p_bitint, n_bitint = all_bitint_prod[idx]
            #println("idx $idx $p_bitint $(int2bitstr(p_bitint, length(p_msps))) $n_bitint $(int2bitstr(n_bitint, length(n_msps)))")
            if (p_bitint in candidates_p) && (n_bitint in candidates_n)
                push!(idxs_subspace, idx)
            end
        end
        if verbose >= 1
            println("idxs_lowest_filling: $(idxs_subspace)")
            for idx in idxs_subspace
                p_bitint, n_bitint = all_bitint_prod[idx]
                println("$(int2bitstr(p_bitint, length(p_msps))) ⊗ $(int2bitstr(n_bitint, length(n_msps)))")
            end
        end
        if only_2p2h
            return idxs_subspace
        end
        if length(idxs_subspace) < maxnum_subspace_basis
            remaining_idxs = setdiff(1:length(all_bitint_prod), idxs_subspace)
            remaining_idxs = sample(remaining_idxs, maxnum_subspace_basis - length(idxs_subspace), replace=false)
            idxs_subspace = vcat(idxs_subspace, remaining_idxs)
        end
        return idxs_subspace
    end
    @error "Invalid sampling scheme: $sampling_scheme. Supported schemes are :uniform, :lowest_filling, and :lowest_filling_2p2h."
end


function prep_Hcm2bmat(Hamil_snt)
    Hcm2bmat = Dict{Tuple{UInt,UInt,UInt,UInt,UInt}, Float64}( )
    for ch in ["pp", "nn", "pn"]
        sps_1 = ch == "nn" ? Hamil_snt.n_sps : Hamil_snt.p_sps
        sps_2 = ch == "pp" ? Hamil_snt.p_sps : Hamil_snt.n_sps
        for i = 1:length(sps_1)
            n_i, l_i, j_i, tz_i = sps_1[i].n, sps_1[i].l, sps_1[i].j, sps_1[i].tz
            tkey_i = get_nkey4_shift(n_i, l_i, j_i, tz_i)
            for j = 1:length(sps_2)
                n_j, l_j, j_j, tz_j = sps_2[j].n, sps_2[j].l, sps_2[j].j, sps_2[j].tz
                Tz_bra = tz_i + tz_j
                tkey_j = get_nkey4_shift(n_j, l_j, j_j, tz_j)
                for k = 1:length(sps_1)
                    n_k, l_k, j_k, tz_k = sps_1[k].n, sps_1[k].l, sps_1[k].j, sps_1[k].tz
                    if tz_k != tz_i
                        continue
                    end
                    tkey_k = get_nkey4_shift(n_k, l_k, j_k, tz_k)
                    for l = 1:length(sps_2)
                        n_l, l_l, j_l, tz_l = sps_2[l].n, sps_2[l].l, sps_2[l].j, sps_2[l].tz
                        if (tz_k + tz_l) != Tz_bra
                            continue
                        end
                        tkey_l = get_nkey4_shift(n_l, l_l, j_l, tz_l)
                        Jmin_bra = div(abs(j_i - j_j), 2)
                        Jmax_bra = div(j_i + j_j, 2)
                        Jmin_ket = div(abs(j_k - j_l), 2)
                        Jmax_ket = div(j_k + j_l, 2)
                        Jmin = max(Jmin_bra, Jmin_ket)
                        Jmax = min(Jmax_bra, Jmax_ket)
                        for Jtot = Jmin:Jmax
                            nkey = (tkey_i, tkey_j, tkey_k, tkey_l, UInt(Jtot))
                            val = hcm_2body(n_i, l_i, j_i, tz_i, n_j, l_j, j_j, tz_j,
                                            n_k, l_k, j_k, tz_k, n_l, l_l, j_l, tz_l,
                                            Jtot)
                            Hcm2bmat[nkey] = val                        
                        end
                    end
                end
            end
        end
    end
    return Hcm2bmat

end

function parse_e3max_from_filename(fn_3nf)
    e1max_file = e2max_file = e3max_file = 0
    #ThBME_srg2.0_ramp40-5-36-7-32-9-28-11-24_N3LO_EM500_c1_-0.81_c3_-3.2_c4_5.4_cD_0.7_cE_-0.06_LNL2_650_500_IS_hw20from30_ms1_2_1.me3j.gz 
    # extract ms_XX_YY_ZZ pattern
    pattern = r"ms(\d+)_(\d+)_(\d+)"
    m = match(pattern, fn_3nf)
    if m !== nothing
        e1max_file = parse(Int, m.captures[1])
        e2max_file = parse(Int, m.captures[2])
        e3max_file = parse(Int, m.captures[3])
        println("Extracted from filename: e1max=$e1max_file, e2max=$e2max_file, e3max=$e3max_file")
    else
        println("Could not extract e3max information from filename: $fn_3nf. Using default e*max=0.")
    end
    return e1max_file, e2max_file, e3max_file
end

function read_3NF_into_snt!(Hamil_snt, fn_3nf, e1max, e3max, to)
    println("Reading 3NF from file: $fn_3nf")
    # Try extract model space from the filename
    e1max_file, e2max_file, e3max_file = parse_e3max_from_filename(fn_3nf)
    sps_p = Hamil_snt.p_sps
    sps_n = Hamil_snt.n_sps
    sps_snt = vcat(sps_p, sps_n)

    # Prealloc Wigner coeffs.
    jmax = 2 * e1max + 1
    dcg_spin = prep_dcg_spin()
    d6j_lj = prep_d6j_lj(jmax)
    dWS = dWS2n(Dict{Int64, Float64}(), Dict{Int64, Float64}(),
                dcg_spin, Dict{UInt64, Float64}(), d6j_lj,
                Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}(),
                Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}())

    _Obj_3BME = main_read_me3j(fn_3nf, 
                               e1max, e1max_file, e2max_file, e3max, 
                               e3max_file, sps_snt, dWS, to) 
                
    push!(Hamil_snt.V3b, _Obj_3BME)
    return nothing 
end


"""
Hrank: rank of the Hamiltonian operator (e.g., 2 for two-body interaction) currently not much used, but may be useful for future extensions.
"""
function qsci_main(sntf, target_nuc, parity, Mtot, Hrank, n_eigen;
                   fn_3nf::AbstractString = "",
                   e1max::Int= 0,
                   e3max::Int= 0,
                   verbose::Int=0,
                   ncsm::Bool=false,
                   is_eval_EntropyMeasure::Bool=false,
                   sampling_method::Symbol=:exact,
                   maxnum_subspace_basis::Int=0,
                   sampled_bits::Union{Vector{Int}, Vector{Int128}, Vector{String}, Nothing}=nothing,
                   get_evecs::Bool=true,
                   ret_evecs::Bool=false,
                   is_show::Bool=false,
                   beta_cm = 0.0,
                   is_RDM::Bool=false,
                   show_breakdown_eigenstates::Bool=false,
                   nu_max::Int=10^5,
                   reorder_orbits::Bool=false
                   )
    @assert sampling_method in [:exact, :qsci, :QSCI, :random, :seniority] "Invalid sampling method: $sampling_method"
    if ncsm && beta_cm == 0.0
        @warn "You are specifying β_cm = 0, which is not recommended for NCSM calculations. Consider using a non-zero value for β_cm."
    end
    if beta_cm != 0.0  
        if !get_evecs 
            @warn "When beta_cm != 0.0, you need to set the optional argument `get_evecs=true` to get center-of-mass corrected energies.\nHence, `get_evecs` is set to true. If it is not what you want, please stop here."
            get_evecs = true
        end
    end
    method_to_sample_bitstrings = "exact" # default method
    if sampling_method in [:qsci, :QSCI]
        method_to_sample_bitstrings = "qsci"
    elseif sampling_method in [:random]
        method_to_sample_bitstrings = "random"
    elseif sampling_method == :seniority || (sampling_method == :exact && nu_max < 20)
        method_to_sample_bitstrings = "seniority"
    end
    println("Using sampling method: $method_to_sample_bitstrings")
    to = TimerOutput()

    # Read the interaction file and prepare the single-particle states
    Hamil_snt = read_smsnt_dev(sntf, target_nuc; ncsm=ncsm, reorder_orbits=reorder_orbits)

    # Read 3NF if provided
    if fn_3nf != ""
        read_3NF_into_snt!(Hamil_snt, fn_3nf, e1max, e3max, to)
    end

    Acore = Hamil_snt.cp + Hamil_snt.cn
    corenuc = element[Hamil_snt.cp] * string(Acore)
    nucleus = def_nuc(target_nuc, "", corenuc)
    vZ = nucleus.Z - nucleus.cZ
    vN = nucleus.N - nucleus.cN
    
    ## This, considering partitions, would be more efficient for future calculations...
    # p_ptn, n_ptn, prod_pn_ptn = gen_partition_from_snt(sntf, parity, vZ, vN, target_nuc)

    # prepare the single-particle states in M-scheme
    p_sps = Hamil_snt.p_sps; n_sps = Hamil_snt.n_sps       
    if length(p_sps) != length(n_sps)
        @warn "The number of proton and neutron single-particle states are different: $(length(p_sps)) != $(length(n_sps)). This may lead to unexpected results with the current implementation."
    end

    p_msps = sps2msps(p_sps)
    n_msps = sps2msps(n_sps)
    dict_m2j, dict_j2m = make_m2j_j2m_dicts(p_sps, n_sps)
    Nref_proton, Nref_neutron = get_Nref_from_nucleus(p_msps, n_msps, vZ, vN)

    # Precalulate the Clebsch-Gordan coefficients
    int_shift, dict_CGs = prepare_CGcoeffs(p_msps, n_msps)

    # Prepare the configurations
    # This could be either an exact calculation or sampled bitstrings by e.g. Monte Carlo, Quantum computing, etc.
    all_bitint_prod, mdim = nothing, 0

    ## Note: current implementation firstly prepares all configurations in the model space, and then samples them.
    ## This is obviously not applicable for larger model spaces. We need to change this in the future.
    org_all_bitint_prod = prepare_all_configs_in_modelspace(parity, Mtot, vZ, vN, p_msps, n_msps, Nref_proton, Nref_neutron, verbose, to)
    org_mdim = length(org_all_bitint_prod)
    maxnum_subspace_basis = ifelse(maxnum_subspace_basis == 0, org_mdim, maxnum_subspace_basis)

    idxs_subspace = Int[ ] 
    Mask = Dict{UInt, Float64}() # Mask for the subspace basis
    if method_to_sample_bitstrings == "exact"
        all_bitint_prod = org_all_bitint_prod
        org_all_bitint_prod = nothing
    elseif method_to_sample_bitstrings == "random"
        Random.seed!(1234)
        idxs_subspace = random_sampling_of_configs(Hamil_snt, dict_j2m, vZ, vN, org_all_bitint_prod, maxnum_subspace_basis, verbose, sampling_scheme=:lowest_filling_2p2h, nu_max=nu_max)
        if length(idxs_subspace) > maxnum_subspace_basis # This can happen when the model space is small
            @warn "Probably unexpected truncation is introduced: length(idxs_subspace)=$(length(idxs_subspace)) > maxnum_subspace_basis=$(maxnum_subspace_basis)."
            idxs_subspace = idxs_subspace[1:maxnum_subspace_basis] 
        end
        for (idx_sub, idx) in enumerate(idxs_subspace)
            nkey = hash_2ints(idx, idx_sub)
            Mask[nkey] = 1.0
        end
        all_bitint_prod = [org_all_bitint_prod[i] for i in idxs_subspace]
    elseif method_to_sample_bitstrings == "qsci"
        # First, we prepare lowest filling + 2p2h configurations
        idxs_subspace = random_sampling_of_configs(Hamil_snt, dict_j2m, vZ, vN, org_all_bitint_prod, maxnum_subspace_basis, verbose, sampling_scheme=:lowest_filling_2p2h, only_2p2h=true, nu_max=nu_max)
        preferred_configs = [org_all_bitint_prod[i] for i in idxs_subspace]
        println("# of configs. within HF2p2h: $(length(preferred_configs))")
        # sampled bitstrings from e.g. IBM Qiskit
        @assert sampled_bits !== nothing "For QSCI sampling, you need to provide the sampled bits as a Vector{Int}, Vector{Int128} or Vector{String}."
        all_bitint_prod = summarize_bitstrings(sampled_bits, p_msps, n_msps, vZ, vN, maxnum_subspace_basis; Qiskit_ordered=true, 
                                               Nref_proton=Nref_proton, Nref_neutron=Nref_neutron,
                                               postselection=[vZ+vN, parity, Mtot],
                                               verbose=verbose>=1,
                                               preferred_configs=preferred_configs)
        if length(all_bitint_prod) > maxnum_subspace_basis # This can happen when the model space is small
            @warn "Probably unexpected truncation is introduced: length(idxs_subspace)=$(length(all_bitint_prod)) > maxnum_subspace_basis=$(maxnum_subspace_basis)."
            all_bitint_prod = all_bitint_prod[1:maxnum_subspace_basis]
        end
        Mask = get_Mask_qsci(org_all_bitint_prod, all_bitint_prod, maxnum_subspace_basis, verbose)
    elseif method_to_sample_bitstrings == "seniority"
        idxs_subspace = random_sampling_of_configs(Hamil_snt, dict_j2m, vZ, vN, org_all_bitint_prod, maxnum_subspace_basis, verbose, sampling_scheme=:seniority, nu_max=nu_max)
        idxs_subspace = idxs_subspace[1:min(maxnum_subspace_basis, length(idxs_subspace))]
        all_bitint_prod = [org_all_bitint_prod[i] for i in idxs_subspace]
        Mask = get_Mask_qsci(org_all_bitint_prod, all_bitint_prod, maxnum_subspace_basis, verbose)
    else
        @error "Unsupported sampling method: $method_to_sample_bitstrings"
    end
    mdim = length(all_bitint_prod)
    println("subdim $mdim")

    Hoffd = Dict{UInt, Float64}( )
    if sampling_method != :exact
        @timeit to "Full Hamil const." FullHmat = construct_Hmat(Hamil_snt, org_mdim, org_all_bitint_prod, p_msps, n_msps, vZ, vN, dict_m2j, int_shift, dict_CGs, verbose, to; Hrank=Hrank)
        Hoffd = get_off_diagonal(FullHmat, idxs_subspace); FullHmat = nothing
    end
    
    @timeit to "Many-body Hamil const." Hmat = construct_Hmat(Hamil_snt, mdim, all_bitint_prod, p_msps, n_msps, vZ, vN, dict_m2j, int_shift, dict_CGs, verbose, to; Hrank=Hrank)
    @timeit to "JJ const." Jmat = construct_Jmat(mdim, all_bitint_prod, p_msps, n_msps, verbose, to)
    
    HCMmat = nothing
    if beta_cm != 0.0
        dict_Hcm = prep_Hcm2bmat(Hamil_snt)
        @timeit to "Hcm const." begin
            HCMmat = construct_HCMmat(beta_cm, mdim, all_bitint_prod, p_msps, n_msps, vZ, vN,
                                    dict_m2j, int_shift, dict_CGs, dict_Hcm, verbose, to; Hrank=Hrank)
        end

        # debug purpose: test HCM operation
        if false
            vec_test_1 = zeros(Float64, mdim)
            vec_test_2 = zeros(Float64, mdim)
            for idx_to_fill in 1:mdim
                vec_test_1 .= 0.0
                vec_test_2 .= 0.0
                vec_test_1[idx_to_fill] = 1.0 # test for a specific basis
                keysvec = collect(keys(HCMmat))
                partials = zeros(Float64, mdim, Threads.maxthreadid())
                operate_H_on_vec!(vec_test_2, HCMmat, vec_test_1, partials, keysvec, to)
                HCMv = copy(vec_test_2)
                Hcm_test = dot(vec_test_1, vec_test_2) / beta_cm

                println("<Hcm> for e_$(idx_to_fill):", @sprintf("%8.3f", Hcm_test))
                for i = 1:mdim
                    n_bitint_ket = all_bitint_prod[idx_to_fill][2]
                    n_str_ket = int2bitstr(n_bitint_ket, length(n_msps))
                    p_bitint, n_bitint = all_bitint_prod[i]
                    n_str = int2bitstr(n_bitint, length(n_msps))
                    tmp =  HCMv[i]/beta_cm
                    if abs(tmp) < 1e-3
                        continue
                    end
                    println("<$n_str ⊗ $n_str_ket>", @sprintf("%3i", i), @sprintf("%10.6f", HCMv[i]/beta_cm))
                end
            end
        end        

        # Add the center-of-mass Hamiltonian to the many-body Hamiltonian
        add_HCM_to_Hmat!(Hmat, HCMmat)
    end   

    @timeit to "Lanczos" evals, evecs = lanczos(Hmat, mdim, n_eigen, get_evecs, to; itnum=300, tol=1e-9, debug_mode=2)
    n_eigen = min(n_eigen, length(evals))

    if is_eval_EntropyMeasure
        calc_entropy_and_mutual_information(sntf, target_nuc, evecs, all_bitint_prod,
                                            length(p_msps)+length(n_msps),
                                            p_msps, n_msps, verbose)
    end

    if show_breakdown_eigenstates
        breakdown_eigenstates(evecs[:,1], all_bitint_prod, p_msps, n_msps, (vZ+vN)%2!=0, verbose)
    end

    if verbose >= 2
        println("Eigenvector configs& coefficients of the lowest state:")
        show_eigenvector_configs_coeffs(evecs[:, 1], all_bitint_prod, p_msps, n_msps)
    end

    Ecms = zeros(Float64, n_eigen)        
    if beta_cm != 0.0 
        @timeit to "<Ecm>" eval_HCMexpec(beta_cm, evecs, mdim, HCMmat, Ecms, n_eigen, to)
    end
    Energies = evals - Ecms
    ExEnergies = Energies .- Energies[1]
    print_vec("Energies (MeV):", Energies)
    print_vec("Ex. (MeV):     ", ExEnergies)

    JJvals = eval_JJexpec(evecs, mdim, Jmat, Mtot, n_eigen, to)
    print_vec("<J^2>:         ", JJvals)
    Jvals = get_Jvals_fromJJ(JJvals)
    print_vec("<J>?:          ", Jvals)
    
    ## <H^2> calculation will be implemented hopefully soon.
    evars = zeros(Float64, n_eigen)
    if sampling_method != :exact
        @timeit to "H2 eval." evars = eval_H2(evecs, evals, Mask, Hoffd, org_mdim, to)        
        print_vec("Energy variance:", evars; ine=true)
    end

    if is_RDM
        verbose = 1
        @timeit to "RDM eval." eval_RDM(evecs, all_bitint_prod, p_msps, n_msps, verbose)
    end

    if !ret_evecs
        evecs = zeros(Float64, 1, 1)
    end
    if is_show
        show(to); 
    end
    println("")

    return Result_QSCI(evals, evecs, Jvals, evars, mdim)
end

function show_eigenvector_configs_coeffs(evecs, all_bitint_prod, p_msps, n_msps)
    for i in eachindex(all_bitint_prod)
        p_bitint, n_bitint = all_bitint_prod[i]
        coeff = evecs[i, 1]
        if abs(coeff) > 1e-6
            if p_bitint == 0
                println("... ⊗ $(int2bitstr(n_bitint, length(n_msps))) coeff = $(@sprintf("%10.6f", coeff))")
            elseif n_bitint == 0
                println("$(int2bitstr(p_bitint, length(p_msps))) ⊗ ... coeff = $(@sprintf("%10.6f", coeff))")
            else
                println("$(int2bitstr(p_bitint, length(p_msps))) ⊗ $(int2bitstr(n_bitint, length(n_msps))) coeff = $(@sprintf("%10.6f", coeff))")
            end 
        end
    end
    return nothing
end

function get_Mask_qsci(org_all_bitint_prod, all_bitint_prod, maxnum_subspace_basis, verbose)
    Mask = Dict{UInt, Float64}()
    if verbose >= 1
        println("Mask for QSCI sampling:")
    end
    for (idx_sub, bitint_prod) in enumerate(all_bitint_prod)
        p_bitint, n_bitint = bitint_prod
        for (idx_org, org_bitint_prod) in enumerate(org_all_bitint_prod)
            if org_bitint_prod == bitint_prod
                nkey = hash_2ints(idx_org, idx_sub)
                Mask[nkey] = 1.0
                if verbose >= 1
                    println("Mask[($(idx_org), $(idx_sub))] = 1.0 for $(int2bitstr(p_bitint, length(p_msps))) ⊗ $(int2bitstr(n_bitint, length(n_msps)))")
                end
                break
            end
        end
    end
    return Mask
end

function breakdown_eigenstates(evec, all_bitint_prod, p_msps, n_msps, odd, verbose)
    @assert length(evec) == length(all_bitint_prod) "The length of evec must match the length of all_bitint_prod."
    prob_sens = zeros(Float64, 10)
    for idx in 1:length(evec)
        if abs(evec[idx]) < 1e-6
            continue
        end
        bitint = all_bitint_prod[idx]
        p_bitint, n_bitint = bitint
        vp, vn = count_seniority(bitint, p_msps, n_msps)
        if verbose >= 2 && length(evec) <= 1000
            println("$(int2bitstr(p_bitint, length(p_msps))) (νp=$vp) ⊗ $(int2bitstr(n_bitint, length(n_msps))) (νn=$vn): $(@sprintf("%10.6f", evec[idx]))")
        end
        idx_sens = div(vp + vn, 2) + 1
        prob_sens[idx_sens] += evec[idx]^2
    end
    labels = ifelse(odd, "g.s. sens: v=1, 3,...", "g.s. sens: v=0, 2,...")
    print_vec(labels, prob_sens)
    return nothing
end

function get_off_diagonal(FullHmat, idx_subspace)
    keysvec = collect(keys(FullHmat))
    partials = [Dict{UInt, Float64}() for _ in 1:Threads.maxthreadid()]
    @threads for idx in 1:length(keysvec)
        tkey = keysvec[idx]
        idx_bra, idx_ket = unhash_2ints(tkey)
        if ( (idx_bra in idx_subspace) && (idx_ket in idx_subspace) ) || ( !(idx_bra in idx_subspace) && !(idx_ket in idx_subspace) )
            continue
        end
        tid = Threads.threadid()
        partials[tid][tkey] = FullHmat[tkey]
    end
    Hoffdiag = Dict{UInt, Float64}()
    for tid in 1:Threads.maxthreadid()
        for (tkey, value) in partials[tid]
            if haskey(Hoffdiag, tkey)
                Hoffdiag[tkey] += value
            else
                Hoffdiag[tkey] = value
            end
        end
    end
    return Hoffdiag
end

function count_seniority(bitint, p_msps, n_msps)
    pbitint, nbitint = bitint
    proton_occ = Int[ ]
    neutron_occ = Int[ ]
    for idx_bit in 1:length(p_msps)
        onehot = Int128(1) << (idx_bit-1)
        if onehot & pbitint == onehot
            push!(proton_occ, idx_bit)
        end
    end
    for idx_bit in 1:length(n_msps)
        onehot = Int128(1) << (idx_bit-1)
        if onehot & nbitint == onehot
            push!(neutron_occ, idx_bit)
        end
    end
    vp = vn = 0
    for pn in 1:2
        occs = ifelse(pn == 1, proton_occ, neutron_occ)
        msps = ifelse(pn == 1, p_msps, n_msps)
        for (idx, idx_p) in enumerate(occs)
            orb = msps[idx_p]
            m = orb.jz
            if m > 0
                continue
            end
            hit = 0
            for idx_2 in idx+1:length(occs)
                idx_p_2 = occs[idx_2]
                orb_2 = msps[idx_p_2]
                if delta_morb_except_m(orb, orb_2) == 0
                    continue
                end
                mp = orb_2.jz
                if mp == -m 
                    hit += 1
                end
            end
            @assert hit == 0 || hit == 1 "hit = $hit"
            if pn == 1
                vp += 2 * (1 - hit)
            else
                vn += 2 * (1 - hit)
            end
        end
    end

    return vp, vn
end