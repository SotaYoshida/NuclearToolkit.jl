module NuclearToolkit

using AssociatedLegendrePolynomials
using Arpack
using Base.Threads
using CodecZlib
using Combinatorics
using DocStringExtensions
using FLoops
using Glob
using HDF5
using LatinHypercubeSampling
using LaTeXStrings
using LinearAlgebra
using Measures
using MKL # for intel machine
using MPI
using Parsers
using Plots
using Printf
using Random
using StatsBase
using Statistics
using SpecialFunctions
using TimerOutputs
using WignerSymbols

include("util_func.jl")
export latex_nuc, int2bitstr, reg, element, rm_comment, rm_nan, print_vec, show_matrix
export get_nkey2_arr, get_nkey2, get_nkey2_u, get_nkey3
export get_nkey3_u, get_nkey3_JPT, get_nkey3_ketJ
export get_nkey4
export get_nkey6, get_nkey6_shift, get_nkey_from_abcdarr
export show_allitems_of_dict
export deltaf

### ChiEFTint.jl  chiEFTint/
include("chiEFTint/struct_const_io.jl")
include("chiEFTint/dict_LECs.jl")
include("chiEFTint/contact.jl")
include("chiEFTint/pionexchange.jl")
include("chiEFTint/angmom_algebra.jl")
include("chiEFTint/eff3nf.jl")
include("chiEFTint/main_chiEFTint.jl")
include("chiEFTint/calibration.jl")
include("chiEFTint/valence.jl")
include("chiEFTint/renorm.jl")
include("chiEFTint/threebodyforce.jl")
include("chiEFTint/threebody_Jacobi.jl")
include("chiEFTint/threebody_JacobiHO.jl")
include("chiEFTint/threebody_lab.jl")
include("chiEFTint/matter.jl")
include("chiEFTint/nn_sampling.jl")
include("chiEFTint/read_me3j.jl")
export print_vec
export make_chiEFTint
export nn_IMSRG_sampling
export test3NF

### NuclData.jl
include("NuclData.jl/amedata.jl")

### HartreeFock.jl
include("hartreefock.jl/def_struct.jl")
export SingleParticleState
include("hartreefock.jl/io_input.jl")
export basedat, def_nuc
export readsnt, readsnt_bin
include("hartreefock.jl/main.jl")
include("hartreefock.jl/hf_mbpt.jl")
include("hartreefock.jl/operator.jl")
export nuclist
export hf_main

### IMSRG.jl
include("IMSRG.jl/imsrg_util.jl")
include("IMSRG.jl/commutator.jl")
include("IMSRG.jl/valencespace.jl")
include("IMSRG.jl/emulator_imsrg.jl")
export imsrg_flow_check
export dmd_main

### ShellModel.jl
include("ShellModel/shellmodel_main.jl")
include("ShellModel/lanczos_methods.jl")
include("ShellModel/transit.jl")
include("ShellModel/input_int_snt.jl")
include("ShellModel/eigenvector_continuation.jl")
include("ShellModel/KSHELL.jl")
include("ShellModel/betadecay.jl")
include("ShellModel/trans_snt_msnt.jl")
include("ShellModel/variational_ansatz.jl")
export main_sm,samplerun_sm # from shellmodel.
export readsmsnt
export construct_msps
export prepEC,solveEC,solveEC_UQ # from eigenvector_continuation.jl
export transit_main # from transit.jl
export read_kshell_summary
export eval_betadecay_from_kshell_log
export main_trans_msnt
include("ShellModel/count_dim.jl")
export count_CIdim
export get_msps_from_jsps
export SingleParticleState_Mscheme
export count_only_num_particle
export get_Nref_from_nucleus

include("ShellModel/partitions.jl")
export gen_partition_from_snt

include("Quantum/qsci.jl")
#export read_bitstr, read_bitint
export qsci_main
end




