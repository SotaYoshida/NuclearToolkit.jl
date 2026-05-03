struct OneBodyRDM
    pp::Matrix{Float64}  # One-body reduced density matrix (1-RDM) for protons
    nn::Matrix{Float64}  # One-body reduced density matrix (1-RDM) for neutrons
end

"""
Evaluate the reduced density matrix (RDM) from the eigenvectors and bitstring products.
"""
function eval_RDM(evecs, 
    all_bitint_prod, 
    p_msps::Vector{SingleParticleState_Mscheme}, 
    n_msps::Vector{SingleParticleState_Mscheme},
    verbose::Int;
    nth::Int=1,
    )
    rdm_1_pp = zeros(Float64, length(p_msps), length(p_msps))
    rdm_1_nn = zeros(Float64, length(n_msps), length(n_msps))
    dim = length(all_bitint_prod)

    target_evec = @view evecs[:, nth]
    @assert dot(target_evec, target_evec) ≈ 1.0 "The eigenvector should be normalized. $(dot(target_evec, target_evec))"
    for idx_bra in 1:dim
        bra_p_bitint, bra_n_bitint = all_bitint_prod[idx_bra]
        coeff_bra = target_evec[idx_bra]
        for idx_ket in 1:dim            
            ket_p_bitint, ket_n_bitint = all_bitint_prod[idx_ket]
            coeff_ket = target_evec[idx_ket]

            # proton part
            if bra_p_bitint *ket_p_bitint != 0
                eval_AdagA!(bra_p_bitint, ket_p_bitint, p_msps,
                            coeff_bra * coeff_ket, rdm_1_pp)
            end

            # neutron part
            if bra_n_bitint * ket_n_bitint != 0
                eval_AdagA!(bra_n_bitint, ket_n_bitint, n_msps,
                            coeff_bra * coeff_ket, rdm_1_nn)
            end

        end
    end

    if verbose > 0
        show_matrix("1-RDM: π", rdm_1_pp)
        show_matrix("1-RDM: ν", rdm_1_nn)
        evals, evecs = eigen(rdm_1_nn)
        print_vec("eigvals(1-RDM: ν):", evals)
        println("evecs(1-RDM: ν):")
        for i = 1:length(evals)
            txt = "@$i: "
            for j = 1:size(evecs, 1)
                if abs(evecs[j, i]) > 1.e-8
                    txt *= " @ $j  $(evecs[j, i])"
                end
            end
            println(txt)
        end
    end

    return OneBodyRDM(rdm_1_pp, rdm_1_nn)
end

"""
Evaluate the matrix element <bra|a^†_p a_q|ket> for the 1-RDM.
1-RDM matrix will be destructively modified.
"""
function eval_AdagA!(bra_bitint, ket_bitint, msps, coeff, rdm)
    # Evaluate the matrix element <bra|a^†_p a_q|ket>
    diff = bra_bitint ⊻ ket_bitint
    if count_ones(diff) > 2
        return nothing
    end

    if count_ones(diff) == 0
        for idx in 1:length(msps)
            if ket_bitint & (2^(idx-1)) != 0 
                rdm[idx, idx] += coeff 
            end
        end
    end

    if count_ones(diff) == 2
        idx_ani = diff & ket_bitint
        idx_bra = diff & bra_bitint
        idx_p = Int(log2(idx_ani))+1
        idx_q = Int(log2(idx_bra))+1
        phase_ani = phase_cre = 1.0

        # The following phase are to be considered by the number of 1s
        # in the bitstring before cre/ani (from the right)
        if idx_p > 1
            phase_ani = count_ones(ket_bitint & (2^(idx_p-1) - 1)) % 2 == 0 ? 1.0 : -1.0
        end
        if idx_q > 1
            phase_cre = count_ones(bra_bitint & (2^(idx_q-1) - 1)) % 2 == 0 ? 1.0 : -1.0
        end
        coeff *= phase_ani * phase_cre

        rdm[idx_p, idx_q] += coeff #* phase_ani * phase_cre
        rdm[idx_q, idx_p] += coeff #* phase_cre * phase_ani

        if ((idx_p, idx_q) == (1, 9)  ) && abs(coeff) > 1e-6
           println(" a1a9: coeff = $(@sprintf("%10.6f", coeff)), ",
                   "bra: $(int2bitstr(bra_bitint, length(msps))), ket: $(int2bitstr(ket_bitint, length(msps))) phase_ket = $phase_ani, phase_bra = $phase_cre")
        end
        if ((idx_p, idx_q) == (2, 10) ) && abs(coeff) > 1e-6
           println("a2a10: coeff = $(@sprintf("%10.6f", coeff)), ",
                   "bra: $(int2bitstr(bra_bitint, length(msps))), ket: $(int2bitstr(ket_bitint, length(msps))) phase_ket = $phase_ani, phase_bra = $phase_cre")
        end
    end

    return nothing
end

# function eval_EntropyMeasure(evecs, all_bitint_prod, 
#     p_msps, n_msps, vZ::Int, vN::Int, verbose::Int;
#     target_states::Vector{Int} = [1])
#     # Evaluate the entanglement entropy measure for the eigenvectors
#     dim = length(all_bitint_prod)
#     n_p = length(p_msps)
#     n_n = length(n_msps)

#     for nth in target_states
#         if nth > size(evecs, 2)
#             println("Warning: target state $nth is out of bounds for the eigenvectors. Skipping.")
#             continue
#         end
#         target_evec = @view evecs[:, nth]
#         @assert dot(target_evec, target_evec) ≈ 1.0 "The eigenvector should be normalized. $(dot(target_evec, target_evec))"
#         rdm = eval_RDM(target_evec, all_bitint_prod, p_msps, n_msps, verbose; nth=nth)

#     end
# end

function _compress_bits(bitint::Int, bit_positions)
    out = 0
    for (i, p) in enumerate(bit_positions)
        out |= ((bitint >> (p - 1)) & 1) << (i - 1)
    end
    return out
end

"""
Reduced density matrix rho_A from a state expanded on `basis_configs`.

Parameters
----------
statevec : 1D AbstractVector{<:Complex}
    Coefficients on the many-body basis `basis_configs`.
basis_configs : AbstractVector{<:Integer}
    Integer-encoded occupation basis (m-scheme single-particle basis).
subsystem : AbstractVector{<:Integer}
    Mode indices (single-particle indices in m-scheme) kept in rho_A.
    Index convention: bit position == mode index (LSB is mode 1).

Returns
-------
rho_A : Matrix{ComplexF64}, size (2^length(subsystem), 2^length(subsystem))
"""
function reduced_density_matrix_from_basis_statevector(statevec, basis_configs, subsystem, p_msps, n_msps)
    subsystem = Tuple(subsystem)
    dimA = 1 << length(subsystem)
    rho_A = zeros(ComplexF64, dimA, dimA)

    submask = 0
    for p in subsystem
        submask |= (1 << (p - 1))
    end

    buckets = Dict{Int, Vector{Tuple{Int, ComplexF64}}}()
    for (cfg, amp) in zip(basis_configs, statevec)
        p_bitint, n_bitint = cfg
        # shift n_bitint to the left by length(p_msps) to combine with p_bitint
        combined_bitint = p_bitint | (n_bitint << length(p_msps))
        a = _compress_bits(combined_bitint, subsystem)  # subsystem index (0-based)
        env_key = combined_bitint & (~submask)          # traced-out part (uncompressed key)
        push!(get!(buckets, env_key, Vector{Tuple{Int, ComplexF64}}()), (a, ComplexF64(amp)))
    end

    for entries in values(buckets)
        for (a1, c1) in entries
            for (a2, c2) in entries
                rho_A[a1 + 1, a2 + 1] += c1 * conj(c2)
            end
        end
    end

    tr = real(sum(diag(rho_A)))
    if tr > 0
        rho_A ./= tr
    end
    return rho_A
end

function von_neumann_entropy(rho; log_base=2, eps=1e-12)
    rho = 0.5 * (rho + rho')
    evals = eigvals(Hermitian(rho))
    evals = clamp.(real.(evals), 0.0, 1.0)
    nz = evals[evals .> eps]
    if isempty(nz)
        return 0.0
    end
    if log_base == 2
        return -sum(nz .* log2.(nz))
    end
    return -sum(nz .* log.(nz)) / log(log_base)
end

"""
Compute single-mode entanglement entropy and pairwise mutual information
from eigenvectors in m-scheme basis.

Parameters
----------
evecs : AbstractMatrix or AbstractVector
    Eigenvectors, size (dim_basis, n_states) or (dim_basis,).
basis_configs : AbstractVector{<:Integer}
    Integer-encoded basis configurations.
n_modes : Integer
    Number of single-particle modes (m-scheme orbitals / qubits).
state_indices : AbstractVector{<:Integer} or nothing
    Which eigenvector columns to evaluate. `nothing` => all columns.
log_base : Real
    Log base for entropy (2 gives bits).

Returns
-------
results : Dict{Int, Dict{String, Any}}
    results[k]["S1"] -> Vector{Float64} size (n_modes,)
    results[k]["MI"] -> Matrix{Float64} size (n_modes, n_modes)
    where k is eigenvector column index.
"""
function calc_entropy_and_mutual_information(
    sntf,
    target_nuc, 
    evecs, 
    basis_configs, 
    n_modes, 
    p_msps,
    n_msps, 
    verbose; 
    state_indices=[4],
    log_base=2
    )
    vecs = ndims(evecs) == 1 ? reshape(evecs, :, 1) : evecs
    if size(vecs, 1) != length(basis_configs)
        throw(ArgumentError("evecs first dimension must match length(basis_configs). Got $(size(vecs, 1)), expected $(length(basis_configs))."))
    end

    if state_indices === nothing
        state_indices = axes(vecs, 2)
    end

    results = Dict{Int, Dict{String, Any}}()
    for k in state_indices
        psi = ComplexF64.(vecs[:, k])
        norm_psi = norm(psi)
        if norm_psi == 0
            throw(ArgumentError("State vector at column $k has zero norm."))
        end
        psi ./= norm_psi

        S1 = zeros(Float64, n_modes)
        for i in 1:n_modes
            rho_i = reduced_density_matrix_from_basis_statevector(psi, basis_configs, [i], p_msps, n_msps)
            S1[i] = von_neumann_entropy(rho_i; log_base=log_base)
        end

        MI = zeros(Float64, n_modes, n_modes)
        for i in 1:n_modes
            for j in (i + 1):n_modes
                rho_ij = reduced_density_matrix_from_basis_statevector(psi, basis_configs, [i, j], p_msps, n_msps)
                Sij = von_neumann_entropy(rho_ij; log_base=log_base)
                mij = S1[i] + S1[j] - Sij
                MI[i, j] = MI[j, i] = max(0.0, mij)
            end
        end
        results[Int(k)] = Dict("S1" => S1, "MI" => MI)
    end

    if verbose > 0
        for (k, res) in results
            println("State index: $k")
            print_vec("Single-mode entropies S1:", res["S1"])
            println("Pairwise mutual information MI:")
            show_matrix("", res["MI"])
        end
    end
    plot_MI_matrix(sntf, target_nuc, results[state_indices[1]]["MI"]; n_p=length(p_msps))
    return results
end


"""

Feb. 20: I will add LaTeX string for nuc in plot...
"""
function plot_MI_matrix(sntf, nuc, MI; n_p=nothing, vmax=0.3)  

    MI = Array{Float64}(MI)
    if ndims(MI) != 2 || size(MI, 1) != size(MI, 2)
        throw(ArgumentError("MI must be a square matrix."))
    end

    n_modes = size(MI, 1)
    n_p = isnothing(n_p) ? (n_modes ÷ 2) : n_p
    n_n = n_modes - n_p
    if n_p < 0 || n_p > n_modes
        throw(ArgumentError("n_p must be in [0, n_modes]."))
    end

    pp = fill(NaN, n_modes, n_modes)
    nn = fill(NaN, n_modes, n_modes)
    pn = fill(NaN, n_modes, n_modes)

    if n_p > 0
        pp[1:n_p, 1:n_p] = MI[1:n_p, 1:n_p]
        pn[1:n_p, (n_p + 1):n_modes] = MI[1:n_p, (n_p + 1):n_modes]
        pn[(n_p + 1):n_modes, 1:n_p] = MI[(n_p + 1):n_modes, 1:n_p]
    end
    if n_n > 0
        nn[(n_p + 1):n_modes, (n_p + 1):n_modes] = MI[(n_p + 1):n_modes, (n_p + 1):n_modes]
    end

    fig = CairoMakie.Figure(size=(800, 600))
    ax = CairoMakie.Axis(
        fig[1, 1];
        xlabel="SPS Index",
        ylabel="SPS Index",
    )

    hm_pp = CairoMakie.heatmap!(ax, pp; colormap=:Reds, colorrange=(0, vmax))
    hm_nn = CairoMakie.heatmap!(ax, nn; colormap=:Blues, colorrange=(0, vmax))
    hm_pn = CairoMakie.heatmap!(ax, pn; colormap=:Purples, colorrange=(0, vmax))

    boundary = n_p + 0.5
    CairoMakie.lines!(ax, [boundary, boundary], [0.5, n_modes + 0.5]; color=:black, linewidth=1)
    CairoMakie.lines!(ax, [0.5, n_modes + 0.5], [boundary, boundary]; color=:black, linewidth=1)

    # Make LaTeX string for the nucleus
    cnuc = latex_nuc(nuc)
    CairoMakie.text!(ax, 0.03, 0.98; text=cnuc, color=:black, 
                    fontsize=36, align=(:left, :top), space=:relative)

    CairoMakie.Colorbar(fig[1, 2], hm_pp; label="MI pp")
    CairoMakie.Colorbar(fig[1, 3], hm_nn; label="MI nn")
    CairoMakie.Colorbar(fig[1, 4], hm_pn; label="MI pn")

    # xticksをしてい
    ax.xticks = (1:4:n_modes, string.(1:4:n_modes))
    ax.yticks = (1:4:n_modes, string.(1:4:n_modes))

    csnt = String(split(sntf, "/") |> last)
    csnt = split(csnt, ".") |> first
    fn = "MI_matrix_$(nuc)_$(csnt).pdf"
    save(fn, fig)
    return nothing
end

