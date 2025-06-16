struct Partition
    occ_jj::Tuple{Vararg{Int}}
    parity::Int
    Nexc::Int
end

struct ProdPNPartition
    p_ptn::Partition
    n_ptn::Partition
end

function get_single_Partition(sps::Vector{Vector{Int}}, num_particles::Int, Nmax::Int)
    nbit = length(sps)  
    max_occ = [ sps_jj[3] + 1 for sps_jj in sps]
    Nexc_ref = 0
    occ_n = 0
    for sps_jj in sps
        n, l, j, tz = sps_jj
        e = 2 * n + l
        for m = -j:2:j
            if occ_n == num_particles
                break
            end
            occ_n += 1
            Nexc_ref += e
        end
    end
    #println("Nexc_ref = $Nexc_ref")

    ls = [ sps[i][2] for i in 1:nbit ]
    ns = [ sps[i][1] for i in 1:nbit ]
    ranges = map(x -> 0:x, max_occ)
    partitions = Partition[ ]
    for occ in Iterators.product(ranges...)
        nocc = sum(occ)
        if nocc != num_particles
            continue 
        end
        lsum = sum(occ[i] * ls[i] for i in 1:nbit)
        parity = lsum % 2 == 0 ? 1 : -1
        Nexc = sum( occ[i]*(2*ns[i]+ls[i]) for i in 1:nbit) - Nexc_ref
        if Nexc > Nmax
            continue
        end
        push!(partitions, Partition(occ, parity, Nexc))
    end
    if isempty(partitions)
        push!(partitions, Partition(Tuple(zeros(Int, nbit)), 1, 0))
    end
    return partitions
end

function show_ptn(ptns::Vector{Partition}, title, verbose=1)
    if verbose > 0
        println("=== $title  ===")
        println("# = $(length(ptns))")
    end
    if verbose > 1
        for ptn in ptns
            println("Occupation: $(ptn.occ_jj), Parity: $(@sprintf("%3i", ptn.parity))")
        end
    end
end

function show_prod_ptn(prod_pn_ptn::Vector{Tuple{Int,Int}}, title, verbose=1)
    if verbose > 0
        println("=== $title ===")
        println("# = $(length(prod_pn_ptn))")
    end
    if verbose > 1
        for idx_pn = 1:length(prod_pn_ptn)
            idx_p, idx_n = prod_pn_ptn[idx_pn]
            println("idx_pn $(@sprintf("%5i", idx_pn)) π $(@sprintf("%3i", idx_p))    ν $(@sprintf("%3i", idx_n))")
        end
    end
end

function get_prod_pn_ptn(p_ptn::Vector{Partition}, n_ptn::Vector{Partition}, parity, Nmax)
    prod_pn_ptn = Tuple{Int, Int}[]
    for (idx_p, p) in enumerate(p_ptn)
        for (idx_n, n) in enumerate(n_ptn)
            if p.parity * n.parity == parity && p.Nexc + n.Nexc <= Nmax
                push!(prod_pn_ptn, (idx_p, idx_n))
            end
        end
    end
    return prod_pn_ptn
end

function make_partitions(p_sps::Vector{Vector{Int}}, n_sps::Vector{Vector{Int}}, proton_number, neutron_number, parity, Nmax)
    """
    p_sps is a vector of jj states for protons consisting of... n, l, j, tz
    Now we have to generate all possible partitions allowed by the given proton_number
    A naive approach is to generate all possible combinations of jj states within N bits
    where N is the number of jj states.
    Among them, the specified number of jj states must be occupied.
    For each such combination, we can count the occupations for odd-number l orbitals,
    giving us the parity information for the partition.
    """
    nbit = length(p_sps)

    p_ptn = get_single_Partition(p_sps, proton_number, Nmax)
    n_ptn = get_single_Partition(n_sps, neutron_number, Nmax)

    show_ptn(p_ptn, "Proton Partitions")
    show_ptn(n_ptn, "Neutron Partitions")

    prod_pn_ptn = get_prod_pn_ptn(p_ptn, n_ptn, parity, Nmax)
    show_prod_ptn(prod_pn_ptn, "Product Partitions (p, n) with parity $(@sprintf("%3i", parity))")

    return p_ptn, n_ptn, prod_pn_ptn 
end

function show_sps(sps::Vector{Vector{Int}})
    for idx = 1:length(sps)
        sps_jj = sps[idx]
        n, l, j, tz = sps_jj
        pn = ifelse(tz==-1, "π", "ν")
        alph_l = chara_l[l+1]
        print("$(pn)$(n)$(alph_l)$(j)/2" * ifelse(idx!=length(sps), ", ", "\n"))
    end
end

function eval_ptn_distance(ptn::Partition, ptn_::Partition)
    occ = ptn.occ_jj
    occ_ = ptn_.occ_jj
    dist = 0
    for i in 1:length(occ)
        dist += abs(occ[i] - occ_[i])
    end
    return dist
end

function eval_ptn_connectable(p_ptn::Vector{Partition}, n_ptn::Vector{Partition}, 
                              prod_pn_ptn::Vector{Tuple{Int, Int}})
    connectable = Dict{String, Vector{Tuple{Int, Int}}}()
    connectable["NN"] = Tuple{Int, Int}[]
    connectable["NN+3NF"] = Tuple{Int, Int}[]

    for i in 1:length(prod_pn_ptn)
        idx_p_i, idx_n_i = prod_pn_ptn[i]
        p_i = p_ptn[idx_p_i]
        n_i = n_ptn[idx_n_i]
        for j in i+1:length(prod_pn_ptn)
            idx_p_j, idx_n_j = prod_pn_ptn[j]
            p_j = p_ptn[idx_p_j]
            n_j = n_ptn[idx_n_j]
            distance_p = eval_ptn_distance(p_i, p_j)
            distance_n = eval_ptn_distance(n_i, n_j)
            distance_summ = distance_p + distance_n
            if distance_summ > 6
                continue
            end
            if distance_summ == 6
                push!(connectable["NN+3NF"], (i, j))
            else
                push!(connectable["NN"], (i, j))
                push!(connectable["NN+3NF"], (i, j))
            end
        end
    end
    return connectable
end

"""

Function to plot the partitions, connectable through nuclear interactions.
The sparsity of the connections is influenced by whether 3NF is included or not.
If the (i, j) pair of partitions is connectable, corresponding space is colored (in a square grid).
"""
function plot_ptn(target_nuc, ptn_connectable, p_ptn, n_ptn, prod_pn_ptn, with3NF::Bool)
    figsize = (300, 300); channels = ["NN"]
    if with3NF
        figsize = (600, 300)
        push!(channels, "NN+3NF") 
    end
    f = Figure(size=figsize, backgroundcolor = :transparent)
    for (nth, channel) in enumerate(channels)
        ax = Axis(f[1, nth], 
                  backgroundcolor = :transparent, aspect = 1)
        xlims!(ax, 1, length(prod_pn_ptn))
        ylims!(ax, 1, length(prod_pn_ptn))
        hidedecorations!(ax) 
        mat = zeros(Float64, length(prod_pn_ptn), length(prod_pn_ptn))
        coln = 1.0
        num_nonzero = length(prod_pn_ptn)
        for j = 1:length(prod_pn_ptn)
            mat[j, j] = coln            
        end
        if haskey(ptn_connectable, channel) 
            if isempty(ptn_connectable[channel])
                continue
            end
            for (i, j) in ptn_connectable[channel]
                mat[i, j] = coln
                mat[j, i] = coln
                num_nonzero += 2
            end
            heatmap!(ax, mat, colorrange = (0., 1.), rasterize = 5)
        end
        title = latex_nuc(target_nuc) 
        text!(ax, 0.05, 0.98, text = title, space = :relative, align = (:left, :top),
              color = :white, transparency = 0.0, fontsize=16)

        ch_fin = ifelse(with3NF, "$(channel)", "valence NN")
        title = L"\mathrm{%$(ch_fin)}"
        text!(ax, 0.25, 0.980, text = title, space = :relative, align = (:left, :top),
              color = :white, transparency = 0.0, fontsize=12)

        p = @sprintf("%.1f", log10(length(prod_pn_ptn)))
        title = L"Dim. partition = $10^{%$(p)}$"
        text!(ax, 0.25, 0.94, text = title, space = :relative, align = (:left, :top),
              color = :white, transparency = 0.0, fontsize=12)
        
        println("Sparsity of $channel: $(@sprintf("%9.2e", (num_nonzero / (length(prod_pn_ptn)^2))))")
    end    
    save("Partition_connectable_$(target_nuc)"*ifelse(with3NF, "_w3NF", "")*".pdf", f)    
    return nothing
end

function gen_partition_from_snt(fn::String, parity::Int, proton_number::Int, neutron_number::Int,
                                target_nuc::String;
                                Nmax::Int = 10^6, with3NF::Bool = false,
                                is_plot::Bool = false, is_show::Bool = false)    
    @assert typeof(parity) == Int && (parity == 1 || parity == -1) "Parity must be either 1 or -1"
    to = TimerOutput()
    # Load the interaction file

    SMobj = readsmsnt(fn, target_nuc)
    p_sps = SMobj.p_sps
    n_sps = SMobj.n_sps

    println("=== Model Space ===")
    show_sps(p_sps)
    show_sps(n_sps)
    p_ptn, n_ptn, prod_pn_ptn = make_partitions(p_sps, n_sps, proton_number, neutron_number, parity, Nmax)
    if is_plot
        @timeit to "connectable"  ptn_connectable = eval_ptn_connectable(p_ptn, n_ptn, prod_pn_ptn)
        @timeit to "plot_ptn" plot_ptn(target_nuc, ptn_connectable, p_ptn, n_ptn, prod_pn_ptn, with3NF)
    end
    if is_show && is_plot
        show(to); println("")
    end
    return p_ptn, n_ptn, prod_pn_ptn
end

"""
Read a bit string for a single species of nucleons (protons or neutrons),
and evaluate Nocc, parity, 
"""
function read_bitstr_single(bitstr::String, m_sps::SingleParticleState_Mscheme;
                            Nexc_ref=0, Qiskit_ordered::Bool=false)
    bitstr = ifelse(Qiskit_ordered, reverse(bitstr), bitstr)
    Nq = length(bitstr)
    @assert Nq == length(sps) "Length of bitstr must match the number of single-particle states"
    Nocc = 0
    Nexc = -Nexc_ref
    parity = 1
    for i = 1:Nq
        if bitstr[i] == '1'
            orbit = m_sps[i]
            e = orbit.e
            l = orbit.l
            parity *=(-1)^l
            j = orbit.j
            jz = orbit.jz
            Nocc += 1
            Nexc += e
        end
    end
end

# function demo_get_ptn()

#     sntf = "interaction_file/usdb.snt"; Zcore = 8; Ncore=8; Acore = Zcore + Ncore
#     proton_number = 4
#     neutron_number = 4
#     Anum = Acore + proton_number + neutron_number
#     parity = 1
#     cnuc = string(Anum) * element[Zcore+proton_number]
#     gen_partition_from_snt(sntf, parity, proton_number, neutron_number, Anum, cnuc)


#     sntf = "interaction_file/sdpf-m.snt"; Acore = 16
#     proton_number = 4
#     neutron_number = 12
#     Anum = Acore + proton_number + neutron_number
#     parity = 1
#     gen_partition_from_snt(sntf, parity, proton_number, neutron_number, Anum, cnuc)


#     sntf = "interaction_file/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw20_emax3_e2max6.kshell.snt"

#     parity = 1
#     Acore = 0
#     proton_number = 3
#     neutron_number = 3
#     Nmax = 10
#     Anum = Acore + proton_number + neutron_number
#     cnuc = element[proton_number] * string(Anum)
#     gen_partition_from_snt(sntf, parity, proton_number, neutron_number, Anum, cnuc; Nmax=Nmax, with3NF=true)

#     parity = 1
#     Acore = 0
#     proton_number = 0
#     neutron_number = 8
#     Nmax = 10
#     Anum = Acore + proton_number + neutron_number
#     cnuc = element[proton_number] * string(Anum)
#     cnuc = "n" * string(Anum)
#     gen_partition_from_snt(sntf, parity, proton_number, neutron_number, Anum, cnuc; Nmax=Nmax, with3NF=true)

#     # sntf = "interaction_file/TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw20_emax3_e2max6.kshell.snt"
#     # ref = "nucl"
#     # corenuc = ""
#     # nuc = def_nuc(cnuc, ref, corenuc)
#     # m_sps,dicts1b,dicts = readsnt(sntf, Anum)
#     # mdim = count_CIdim(nuc, m_sps, Nmax)
#     # println("mdim $mdim")

# end
# demo_get_ptn()