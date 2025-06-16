struct SingleParticleState_Mscheme
    e::Int64
    n::Int64
    l::Int64
    j::Int64
    jz::Int64
    tz::Int64
end

function get_msps_from_jsps(sps::Array{SingleParticleState}; concatenate=false)
    msps_proton = SingleParticleState_Mscheme[]
    msps_neutron = SingleParticleState_Mscheme[]
    for i in 1:length(sps)
        n = sps[i].n
        l = sps[i].l
        j = sps[i].j
        tz = sps[i].tz
        e = 2*n + l
        for jz in -j:2:j
            target = ifelse(sps[i].tz == -1, msps_proton, msps_neutron)
            push!(target, SingleParticleState_Mscheme(e,n,l,j,jz,tz))
        end
    end
    if concatenate
        msps = vcat(msps_proton, msps_neutron)
        return msps
    else
        return msps_proton, msps_neutron
    end
end

"""
Counting parity, Nmax, and Mtot for a given configuration. 
Nexc is the quanta for excitations across major shell gaps, and Mtot is the total angular momentum projection.
It should be noted that (m)sps are in ascending order in terms of energy, while the bit string (Integer) is in descending order,
i.e. rightmost bit corresponds to the lowest single particle state.
"""
function count_parity_Nmax_Mtot(bint::Int128, msps, Nref)
    parity = 1
    Nexc = 0
    Mtot = 0
    #println("bint $bint\nbit $(int2bitstr(bint, length(msps)))")
    for idx in 1:length(msps)
        if (bint & (Int128(1) << (idx-1))) != 0
            parity *= (-1)^(msps[idx].l)
            Nexc += msps[idx].e
            Mtot += msps[idx].jz
            #println("hot idx $j Nexc $Nexc {nlj,jz,e} $(msps[idx].e) $(msps[idx].l) $(msps[idx].j) $(msps[idx].jz) $(msps[idx].e)")
        end
    end
    Nexc = Nexc - Nref
    return parity, Nexc, Mtot
end


"""
Counting the sum of e for naive filling configurations, which is needed to eval excitation quanta.
"""
function get_Nref_from_nucleus(p_msps, n_msps, Z, N)
    Nref_proton = sum( [ msps.e for msps in p_msps[1:Z]] )
    Nref_neutron= sum( [ msps.e for msps in n_msps[1:N]] )
    return Nref_proton, Nref_neutron
end

function extract_ji_jj(bint, msps)
    jjs = Int[ ]
    for j in 1:length(msps)
        if (bint & (1 << (j-1))) != 0
            idx = length(msps) + 1 - j
            push!(jjs, msps[idx].j)
        end
    end
    @assert length(jjs) == 2 "Something wrong at extract_ji_jj!! bint = $bint jjs = $jjs"
    return jjs
end

struct config_struct
    Mtot::Int
    parity::Int
    Nexc::Int
    configs::Vector{Int128}
    count::Vector{Int128}
end

"""
To get the location of most significant bit 
"""
function msb_position(x::Integer)
    x == 0 && return 0 
    return sizeof(x) * 8 - leading_zeros(x)
end

function get_idx_hole(e_target, msps)
    te = 0
    idx_r = 1
    while te < e_target
        te = msps[idx_r].e
        idx_r += 1
    end
    return idx_r
end

"""
int_init corresponds to initial integer obtained by naive filling from right most bit.
"""
function get_maxint(int_init, msps, Nmax, verbose=0) ::Int128
    msb = msb_position(int_init)    
    Norb = length(msps)
 
    Nexc = 0
    idxs_annihi = Int[ ]
    idxs_create = Int[ ]
    
    e_Fermi = e_target = msps[msb].e
    # to get rightmost bit of Fermi surface (in terms of emax)
    idx_r = get_idx_hole(e_target, msps)

    for idx_set_1 in Norb:-1:msb+1
        if int_init & (Int128(1) << (idx_r-1)) == 0
            # if idx_r is vacant, idx_r should be determined from lower emax orbitals
            idx_r = get_idx_hole(e_target, msps)
        end
        e_part = msps[idx_set_1].e
        e_hole = msps[idx_r].e
        if e_part - e_hole > Nmax
            continue
        end
        push!(idxs_annihi, idx_r)
        idx_r += 1
        push!(idxs_create, idx_set_1)
        Nexc += e_part - e_hole
        if Nexc == Nmax
            break
        end
    end
    if verbose >= 1
        println("idxs_annihi $idxs_annihi idxs_create $idxs_create")
    end
    x = int_init
    for d in idxs_annihi
        a = Int128(2)^(d-1)        
        x = x ⊻ a
        if verbose >= 1
            println(" ax  ", int2bitstr(x, length(msps)), " x $x")
        end
    end
    for d in idxs_create
        c = Int128(2)^(d-1)
        x = x ⊻ c
    end
    if verbose >= 1
        println("ccx  ", int2bitstr(x, length(msps)), " x $x")
    end
    return x 
end

function hash_3int(Mtot, parity, Nmax)::UInt
    Mtot = Mtot + 1000
    parity = parity + 10
    Nmax = Nmax + 10
    return UInt(Mtot) << 32 | UInt(parity) << 16 | UInt(Nmax)    
end

function unhash_Uint_3int(hash::UInt)::Tuple{Int,Int,Int}
    Mtot = Int(hash >> 32) - 1000
    parity = Int((hash >> 16) & 0xFFFF) - 10
    Nmax = Int(hash & 0xFFFF) - 10
    return (Mtot, parity, Nmax)
end

"""
If xinit(xmaxint) is set to naive filling (largest integer with Nmax+1),
one should wait very long time to get the result.
The following function generates xinit and xmaxint, which is to be distributed to multiple threads.
"""
function generate_xinit_xmax(num_particle, msps, upperbound)
    num_cycles = num_particle
    xinit = Int128(1) << (num_particle) - 1
    max_shift = length(msps) - num_particle 
    x_list = Vector{Int128}[ ]
    basenum = div(num_particle, num_cycles)
    # devide num_particle into num_cycles and shift only them
    nums_bit_shift = [ 0 for _ in 1:num_cycles]
    for i in 1:num_cycles
        nums_bit_shift[i] = i * basenum
    end
    if num_particle % num_cycles != 0
        nums_bit_shift[num_cycles] += num_particle % num_cycles
    end
    ## println("x0:  ", int2bitstr(xinit, length(msps)))
    ## println("xm:  ", int2bitstr(upperbound, length(msps)))
    #println("nums_bit_shift $nums_bit_shift")
    @assert nums_bit_shift[end] == num_particle "nums_bit_shift[end] $(nums_bit_shift[end]) != num_particle"
    for n in 0:max_shift
        x_low = x_R = xinit << n
        if x_R >= upperbound
            break
        end
        msb = msb_position(x_R)
        mlb = msb - num_particle + 1
        for cycle in 1:num_cycles
            n_bit_involved = 0
            if cycle != 0
                n_bit_involved = nums_bit_shift[cycle]
            end
    
            x_L = x_low << 1
            x_L -= x_L & (Int128(1) << (msb - n_bit_involved))

            x_L += (Int128(1) << (mlb -1)) 
            if cycle == num_cycles
                x_L = x_low << 1
            end
            if x_L > upperbound 
                x_L = upperbound
            end 
            #println("$cycle R: ", int2bitstr(x_R, length(msps)))
            #println("$cycle L: ", int2bitstr(x_L, length(msps)), " mlb $mlb")
            # @assert x_L > x_R "n $n cycle $cycle msb $msb n_bit_involved $n_bit_involved x_L $x_L x_R $x_R"

            push!(x_list, [x_R, x_L])
            x_R = x_L
        end
    end
    return x_list
end

function count_only_num_particle(num_particle, msps, Nref, Nmax_specified=10^5;
                                 store_config=false, verbose=0) 
    naive_config = config_struct[ ]
    if num_particle == 0
        push!(naive_config, config_struct(0, 1, 0, [Int128(0)], [Int128(1)]))
    elseif num_particle == length(msps) # fully occupied   
        x = Int128(1) << (num_particle) - 1
        parity, Nmax, Mtot = count_parity_Nmax_Mtot(x, msps, Nref)   
        push!(naive_config, config_struct(Mtot, parity, Nmax, [Int128(1) << (num_particle) - 1], [Int128(1)]))
    else
        int_init = Int128(1) << (num_particle) - 1
        upper_bound = int_init << (length(msps) - num_particle) #???

        x_list = generate_xinit_xmax(num_particle, msps, upper_bound)
        all_encountered = [ Dict{UInt,Vector{Int128}}() for _ in 1:length(x_list) ]
        all_counts = [ Dict{UInt,Vector{Int128}}() for _ in 1:length(x_list) ]

        #println("upper_bound $upper_bound $(int2bitstr(upper_bound, length(msps)))")
        #println("x_list $x_list $(x_list[end])")
        for idx = 1:length(x_list)
            encountered = all_encountered[idx]
            counts = all_counts[idx]
            x = x_list[idx][1]
            maxint = x_list[idx][2] + ifelse(idx==length(x_list), 1, 0)

            while x < maxint
                parity, Nmax, Mtot = count_parity_Nmax_Mtot(x, msps, Nref)   
                if Nmax <= Nmax_specified 
                    tkey = hash_3int(Mtot, parity, Nmax)
                    if haskey(counts, tkey)
                        if store_config
                            push!(encountered[tkey], x)
                        end
                        counts[tkey] .+= Int128(1)
                    else
                        counts[tkey] = [Int128(1)]
                        if store_config
                            encountered[tkey] = [x]
                        end
                    end
                end
                u = x & (-x)
                v = x + u   
                x = v + (((v ⊻ x) ÷ u) >> 2) 
            end
        end

        # First we collect the keys and values from all dictionaries
        all_keys = UInt[ ]
        for i in 1:length(all_counts)
            for key in keys(all_counts[i])
                if !(key in all_keys)
                    push!(all_keys, key)
                end
            end
        end

        for key in all_keys
            Mtot, parity, Nmax = unhash_Uint_3int(key)
            tmp_counts = Int128[0]
            configs = Int128[ ]
            for i in 1:length(all_counts)
                if haskey(all_counts[i], key)
                    tmp_counts[1] += all_counts[i][key][1]                    
                end
                if store_config && haskey(all_encountered[i], key)
                    append!(configs, all_encountered[i][key])
                end
            end
            push!(naive_config, config_struct(Mtot, parity, Nmax, configs, tmp_counts))
        end
    end
    return naive_config
end

function count_CIdim(nuc::nucleus, sps::Array{SingleParticleState},
                     Nmax_specified::Int64;
                     Mspecified::Int64=-1,
                     truncation_scheme::String="FCI(Nmax)",
                     target_parity::Int64=1,
                     verbose::Int=0,
                     show_time=false
    )
    to = TimerOutput()
    Mspecified = if Mspecified == -1        
        nuc.A % 2
    else
        Mspecified
    end

    if truncation_scheme =="FCI(Nmax)"
        @assert Nmax_specified >= 0 "FCI(Nmax) truncation scheme requires argument Nmax_specified > 0 in count_CIdim function"
    end

    # get m-scheme sps from jj-coupling sps
    p_msps, n_msps = get_msps_from_jsps(sps)
    if verbose >= 1
        println("len(p_msps) $(length(p_msps)) len(n_msps) $(length(n_msps))")
        println("Generating configurations with specified Z/N...")
        println("upper value for p/n ", binomial(big(length(p_msps)), nuc.Z), " ", binomial(big(length(n_msps)), nuc.N))
    end
    Nref_proton, Nref_neutron = get_Nref_from_nucleus(p_msps, n_msps, nuc.Z, nuc.N)

    @timeit to "Zconfig" struct_proton_config  = count_only_num_particle(nuc.Z, p_msps, Nref_proton, Nmax_specified)
    @timeit to "Nconfig" struct_neutron_config = count_only_num_particle(nuc.N, n_msps, Nref_neutron, Nmax_specified)

    Nmax_of_interest = 2:2:Nmax_specified

    #println("size of proton/neutron partitions: $(length(struct_proton_config)) $(length(struct_neutron_config))")
    # for i in 1:length(struct_proton_config)
    #     println("dim_p_$i = $(struct_proton_config[i].count[1]) Mtot_p = $(struct_proton_config[i].Mtot) parity_p = $(struct_proton_config[i].parity) Nexc_p = $(struct_proton_config[i].Nexc)")
    # end    
    # for i in 1:length(struct_neutron_config)
    #     println("dim_n_$i = $(struct_neutron_config[i].count[1]) Mtot_n = $(struct_neutron_config[i].Mtot) parity_n = $(struct_neutron_config[i].parity) Nexc_n = $(struct_neutron_config[i].Nexc)")
    # end

    ## To get possible combinations of proton/neutron configurations and count the dimension
    dim_dict = Dict{String, Vector{Int128}}( )
    dim_dict["Nmax"] = Nmax_of_interest
    dim_dict["dim"] = Int128[0 for _ in 1:length(Nmax_of_interest)]

    if truncation_scheme == "FCI(Nmax)"
        num_of_Nmax = length(Nmax_of_interest)
        dim_p_par = length(struct_proton_config)
        dim_n_par = length(struct_neutron_config)
        block_dims = [ zeros(Int128, dim_p_par, dim_n_par) for _ in 1:num_of_Nmax ]

        @timeit to "count" for idx_Nmax = 1:num_of_Nmax
            Nmax_target = Nmax_of_interest[idx_Nmax]
           
            if verbose >= 1
                println("dim_p_par = $dim_p_par dim_n_par = $dim_n_par")
            end
            for idx_pn_prod in 1:dim_p_par*dim_n_par
                idx_p = div(idx_pn_prod, dim_n_par) + ifelse(idx_pn_prod % dim_n_par == 0, 0, 1)
                idx_n = rem(idx_pn_prod, dim_n_par) + 1

                Mtot_p   = struct_proton_config[idx_p].Mtot
                parity_p = struct_proton_config[idx_p].parity
                Nexc_p   = struct_proton_config[idx_p].Nexc
                dim_p = struct_proton_config[idx_p].count[1]
                Mtot_n   = struct_neutron_config[idx_n].Mtot
                parity_n = struct_neutron_config[idx_n].parity
                Nexc_n   = struct_neutron_config[idx_n].Nexc
                dim_n = struct_neutron_config[idx_n].count[1]
                if Mtot_p + Mtot_n != Mspecified
                    continue
                end
                if parity_p * parity_n != target_parity
                    continue
                end
                if Nexc_p + Nexc_n > Nmax_target
                    continue
                end
                block_dims[idx_Nmax][idx_p, idx_n] = dim_p * dim_n
            end
        end
        for idx_Nmax = 1:num_of_Nmax
            Nmax_target = Nmax_of_interest[idx_Nmax]
            dim = sum(block_dims[idx_Nmax])
            dim_dict["dim"][idx_Nmax] = dim
            println("nuc $(@sprintf("%6s", nuc.cnuc))   A = $(@sprintf("%4i", nuc.A))   Z = $(@sprintf("%4i", nuc.Z))   2*Mtot = $(@sprintf("%3i", Mspecified))",
            "  parity = ", ifelse(target_parity == 1, "+", "-"),
            "  Nmax = $(@sprintf("%4i", Nmax_target))  dim(FCI) = 10^($(@sprintf("%4.1f", log10(dim)))) ($dim)")
        end
    else
        @error "Truncation scheme $truncation_scheme not implemented yet"
    end

    if show_time
        show(to; compact=true, allocations =true); print("\n")
    end
    return dim_dict
end

# function CIdim_plot(emax::Int, Data, nucs, Nmaxmax)
#     plot_linewidth = 2
#     marker_pool = [:circle, :diamond, :utriangle, :dtriangle, :pentagon, :hexagon, :rect, ]
#     f = Figure( size=(500, 300), backgroundcolor = :transparent)
#     ax = Axis(f[1, 1] , xlabel = "Nmax", ylabel = "Dim.", yscale = log10,
#               xticks = 2:2:Nmaxmax,
#               yticks = LogTicks(2:2:18),
#               backgroundcolor = :transparent)
#     xlims!(ax, 1, Nmaxmax+1)
#     ylims!(ax, 1e1, 1e18)
#     for idx = 1:length(nucs)
#         marker = marker_pool[idx]
#         nuc = nucs[idx]
#         lnuc = latex_nuc(nuc)
#         tData = Data[nuc]
#         x = tData["Nmax"]
#         y = tData["dim"]
#         if nuc == nucs[1]
#             obj = scatterlines!(ax, x, y,
#                                 label=lnuc, linewidth=plot_linewidth, 
#                                 marker=marker, markersize=15,
#                                 colormap =:seaborn_dark
#             )
#         else
#             obj = scatterlines!(ax, x, y,
#                 label=lnuc, linewidth=plot_linewidth,
#                 marker=marker, markersize=15,
#             )
#         end
#     end
#     axislegend(ax, merge = true, position = (0.01,0.99), fontsize=35, nbanks=2)
#     save("nuclear_model_space_dimensions_emax$(emax).pdf", f)
# end

# function main_count(; plot_from_jld2 = false, show_time = true, log_level = 0)
#     ref = "nucl"
#     corenuc = ""
#     hw = 20
#     emax = 8
#     emax_calc = 4
#     sntf = homedir() * "/Desktop/ThBMEs/tbme_em500n3lo_srg2.0hw$(hw)emax$(emax).snt.bin"
#     println("emax $emax_calc hw $hw ")

#     nucs = ["H3", "He4", "Li6", "n8", "He8", "C12", "O16"]
#     Nmax_max = 16

#     nucs = ["He4", "Li6", "C12", "O16", "Ne20", "Ca40"]
#     nucs = ["He4", "Li6", "Be8", "C12", "O16", "Ne20" ]
#     Nmax_max = 16

#     Data = Dict{String, Dict{String, Vector{Int128}}}( )
#     for cnuc in nucs
#         nuc = def_nuc(cnuc, ref, corenuc)
#         binfo = basedat(nuc, sntf, hw, emax_calc, ref)
#         sps,dicts1b,dicts = readsnt_bin(sntf, binfo; use_Float64=false, neutron_drop=false)
#         Data[cnuc] = count_CIdim(nuc, sps, Nmax_max, show_time=show_time, verbose=log_level)        
#     end

#     # Write out the Data in jld2 format
#     @save "nuclear_model_space_dimensions_emax$(emax_calc).jld2" Data

#     # plot
#     plot_from_jld2 = !false 
#     if plot_from_jld2
#         @load "nuclear_model_space_dimensions_emax$(emax_calc).jld2" Data
#     end

#     CIdim_plot(emax_calc, Data, nucs, Nmax_max)

# end