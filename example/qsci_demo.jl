using TimerOutputs
using HDF5
using CairoMakie
using ColorSchemes
colors = ColorSchemes.seaborn_dark[1:10]

#include("src/NuclearToolkit.jl")
#using .NuclearToolkit
using NuclearToolkit

function read_Qiskit_bits(file_path::String)
    h5f = h5open(file_path, "r") 
    bits = read(h5f, "bits") 
    close(h5f)
    return bits
end

function plot_qsci_demo(target_nuc, Data)
    fig = Figure( size = (600, 400))
    ax = Axis(fig[1, 1], xlabel = "Subspace Dimension", ylabel = "Energy (MeV)")
    
    Data_Exact = Data["Exact"]
    Data_Random = Data["Random"]
    Data_QSCI = Data["QSCI"]
    mdim = collect(keys(Data_Exact))[1]
    evals_Exact = Data_Exact[mdim]
    n_states = length(evals_Exact)

    subdims_R = sort(collect(keys(Data_Random)))
    subdims_Q = sort(collect(keys(Data_QSCI)))
    y_R = [ [ Data_Random[subdim][i] for subdim in subdims_R] for i in 1:n_states ]
    y_Q = [ [ Data_QSCI[subdim][i] for subdim in subdims_Q] for i in 1:n_states ]
    for i in 1:n_states        
        lines!(ax, [0, mdim], [evals_Exact[i], evals_Exact[i]], 
              label = ifelse(i==1,"Exact",nothing), color = colors[i], linestyle=:dot, transparency=0.5)
        scatter!(ax, mdim, evals_Exact[i], label=ifelse(i==1,"Exact",nothing), marker = :star4, color = colors[i])
    end
    for i in 1:n_states        
        scatterlines!(ax, subdims_R, y_R[i], label = ifelse(i==1,"Random",nothing), marker=:circle, color = colors[i], linestyle = :dash)
        scatterlines!(ax, subdims_Q, y_Q[i], label = ifelse(i==1,"TE+QSCI", nothing),  marker=:utriangle, color = colors[i])
    end
    axislegend(ax, position = :lb, merge = true)
    text!(ax, 0.03, 0.98; text = latex_nuc(target_nuc), 
          color = :black, fontsize = 20, align = (:left, :top),
          space = :relative)
    fn = "figures/qsci_demo_$(target_nuc).pdf"
    if !isdir("figures")
        mkpath("figures")
    end
    save(fn, fig)
end

function qscidemo(; verbose::Int=0, Hrank::Int=2, n_eigen::Int=10, nat_parity::Bool=true)

    # Specify the target system, interaction file, parity, and M (total angular momentum projection)
    target_nuc = "O20"
    target_nuc = "F20"
    target_nuc = "Si28"
    target_nuc = "Mg24"

    #target_nuc = "Ne20"
    sntf = "interaction_file/usdb.snt"

    #target_nuc = "Be8"; sntf = "interaction_file/ckpot.snt"

    if !isfile(sntf)
        error("Interaction file $sntf does not exist. Please provide a valid interaction file.")
    end

    # Parity and Mtot are determined from the target_nuc and nat_parity::Bool
    reg = r"[0-9]+"
    Anum = parse(Int, match(reg, target_nuc).match)
    parity = ifelse(Anum % 2 == 0, 1, -1)
    Mtot = ifelse(Anum % 2 == 0, 0, 1)
    if !nat_parity
        parity *= -1
    end

    # Collecting Exact/Random/QSCI results
    Data = Dict{String, Dict{Int, Vector{Float64}}}()
    Data["Exact"] = Dict{Int, Vector{Float64}}()
    Data["Random"] = Dict{Int, Vector{Float64}}()
    Data["QSCI"] = Dict{Int, Vector{Float64}}()
   
    # Call the main API for Exact calculation   
    Res_Exact = qsci_main(sntf, target_nuc, parity, Mtot, Hrank, n_eigen, verbose; 
                          sampling_method = :exact)
    Data["Exact"][Res_Exact.mdim] = Res_Exact.evals

    # Function to read bitstrings from a file
    fn = "sampled_Qiskit_bitstr/qiskit_bitstrings_$(target_nuc)_T1followedbypn.h5"
    #fn = "sampled_Qiskit_bitstr/qiskit_bitstrings_$(target_nuc).h5"
    if !isfile(fn)
        error("File $fn does not exist. Please run the Qiskit sampling script first.")
    end
    sampled_bits = read_Qiskit_bits(fn)

    sub_dims = [ div(Res_Exact.mdim, 10) *i for i in 1:9 ]
    for maxnum_subspace_basis in sub_dims
        # Random CI
        Res_RCI = qsci_main(sntf, target_nuc, parity, Mtot, Hrank, n_eigen, verbose;
                            sampling_method = :random, maxnum_subspace_basis = maxnum_subspace_basis)
        Data["Random"][Res_RCI.mdim] = Res_RCI.evals

        # call the main API for QSCI
        Res_QSCI = qsci_main(sntf, target_nuc, parity, Mtot, Hrank, n_eigen, verbose;
                            sampling_method = :QSCI, maxnum_subspace_basis = maxnum_subspace_basis, 
                            sampled_bits=sampled_bits)
        Data["QSCI"][Res_QSCI.mdim] = Res_QSCI.evals
    end
    
    plot_qsci_demo(target_nuc, Data)
    
end

qscidemo()