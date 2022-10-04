"""
struct `nuclei`
# Fields
- `Z::Int64` proton number of the reference nucleus
- `N::Int64` neutron number of the ref.
- `A::Int64` mass number of the ref.
- `el::String` element (e.g., "He")
- `cnuc::String` string element nameA (e.g., "He8")
- `cZ::Int64` proton number of core nucleus 
- `cN::Int64` neutron number of core
- `corenuc::String` core nucleus (e.g., "He4")
"""  
struct nuclei
    Z::Int64
    N::Int64
    A::Int64
    Aref::Int64
    el::String
    cnuc::String
    cZ::Int64
    cN::Int64
    corenuc::String
end 

"""
struct `chan1b` 
# Fields
- `chs1b::Vector{Dict{Int64,Vector{Int64}}}` dict of single particle states with non-zero contribution (having same l,j) [dict for proton sps, dict for neutron sps]
- `chs1b_redundant::Vector{Dict{Int64,Vector{Int64}}}` redundant version of chs1b (with i>j case)
- `snt2ms::Dict{Int64,Int64}` map from snt idx to modelspace(ms) idx
- `ms2snt::Dict{Int64,Int64}` map from ms idx to snt idx
"""  
struct chan1b 
    chs1b::Vector{Dict{Int64,Vector{Int64}}}
    chs1b_redundant::Vector{Dict{Int64,Vector{Int64}}}
    snt2ms::Dict{Int64,Int64}
    ms2snt::Dict{Int64,Int64}
end

"""
struct `Dict1b` 
# Fields
- `snt2ms::Dict{Int64,Int64}` map from snt idx to modelspace(ms) idx
- `ms2snt::Dict{Int64,Int64}` map from ms idx to snt idxdef_struct.jl
"""
struct Dict1b 
    snt2ms::Dict{Int64,Int64}
    ms2snt::Dict{Int64,Int64}
end

"""
struct `VdictCh`
# Fields
- `Vch::Int64` two-body channel (specified by JPT)
- `Vdict::Dict{Int64,Int64}` dict to get idx from ket, which is used in only `vPandya` function for HFMBPT
"""
struct VdictCh
    Vch::Int64
    Vdict::Dict{Int64,Int64}
end

"""
struct `chan2b`
referred to as "tbc" (two-body channel) in some functions
# Fields
- `Tz::Int64` total tz, -2(pp),0(pn),2(n)
- `prty::Int64` parity
- `J::Int64` total J
- `kets::Vector{Vector{Int64}}` vector of ket (e.g. [1,1], [1,3],...)
"""
struct chan2b
    Tz::Int64
    prty::Int64
    J::Int64
    kets::Vector{Vector{Int64}}
end

"""
struct `Chan2bD`
# Fields
- `Chan2b::Vector{chan2b}` array of chan2b (ch=1,...,nchan)
- `dict_ch_JPT::Dict{Vector{Int64},VdictCh}` dict to get VdictCh by given key `[J,prty,T]`
- `dict_ch_idx_from_ket::Vector{Vector{Dict{Vector{Int64},Vector{Vector{Int64}}}}}` dict to get [ch,idx], having array structure [pnrank(=1/2/3)][J+1], `key`=ket
- `dict_idx_from_chket::Vector{Dict{Vector{Int64},Int64}}` dict to get idx from ket, having array structure [ch]
"""
struct Chan2bD
    Chan2b::Vector{chan2b}
    dict_ch_JPT::Dict{Vector{Int64},VdictCh}
    dict_ch_idx_from_ket::Vector{Vector{Dict{Vector{Int64},Vector{Vector{Int64}}}}}
    dict_idx_from_chket::Vector{Dict{Vector{Int64},Int64}}   
end

"""
struct `single_util122`, used to make operation related commutator122 matrix manipulation
"""
struct single_util122
    ind1_ia::Vector{Int64}
    ind1_ja::Vector{Int64}
    ind2s::Vector{Int64}
    factor_ia::Matrix{Float64}
    factor_ja::Matrix{Float64}
end

"""
struct `PandyaObject`, used for Pandya transformation (especially in `comm222ph_ss!`)
"""
struct PandyaObject
    numbers::Vector{Vector{Int64}}
    numbers_addinv::Vector{Vector{Int64}}
    Chan2b::Vector{chan2b}
    phkets::Vector{Vector{Int64}}
    dict_ich_idx_from_ketcc::Vector{Dict{Int64,Int64}}
    XYbars::Vector{Vector{Matrix{Float64}}}
    Zbars::Vector{Matrix{Float64}}
    PhaseMats::Vector{Matrix{Float64}}
    tMat::Vector{Matrix{Float64}}
    dict_ch2ich::Dict{Int64,Int64}
    keys6j::Vector{Vector{Int64}}
    util122::Vector{Vector{single_util122}}
    Mats_hh::Vector{Matrix{Float64}}
    Mats_pp::Vector{Matrix{Float64}}
    Mats_ph::Vector{Matrix{Float64}}
end


mutable struct valDictMonopole
    monopole::Vector{Float64}
    vals::Vector{Vector{Int64}}
end
"""
struct `dictTBMEs` contains dictionaries for TBME/monopole
# Fields  
- `dictTBMEs::Vector{Dict{Vector{Int64},Float64}}` one can get pp/pn/nn dict by `dictTBMEs[pnrank]` (`pnrank=1,2,3`)
- `dictMonopole::Vector{Dict{Vector{Int64},valDictMonopole}}` one can get monopole component of two-body interaction by `dictMonopole[pnrank][key]`, `key` to be ket array like `[1,1]`
"""
struct dictSnt
    dictTBMEs::Vector{Dict{Vector{Int64},Vector{Float64}}}
    dictMonopole::Vector{Dict{Vector{Int64},valDictMonopole}}
end


"""
struct `basedat` contains base infomation of the calculation
# Fields
- `nuc::nuclei` information of target/core nucleus 
- `sntf::String` filename/path to input interaction
- `hw::Int64` hbar omega parameter used to define single particle states
- `emax::Int64` emax truncation for the entire calculations
- `ref::String` to specify ref="core" or ref="nucl"
"""
struct basedat
    nuc::nuclei
    sntf::String
    hw::Int64
    emax::Int64
    ref::String 
end

"""
struct `hfdata`, used to calculate multiple nuclei in a single runscript
# Fields
- `nuc::nuclei` information of target/core nucleus
- `data::Vector{Vector{Float64}}` will be experimental data from AME2020 (if available)
- `datatype::Vector{String}` supposed to be ["E"] for now
"""
mutable struct hfdata
    nuc::nuclei
    data::Vector{Vector{Float64}}
    datatype::Vector{String}
end

"""
mutable struct `space_channel`;
dictionaries to get the two-body channels that have kets (specified by pp,ph, etc.)
# Fields
- `pp::Dict{Int64,Vector{Int64}}` particle-particle
- `ph::Dict{Int64,Vector{Int64}}` particle-hole 
- `hh::Dict{Int64,Vector{Int64}}` hole-hole
- `cc::Dict{Int64,Vector{Int64}}` core-core
- `vc::Dict{Int64,Vector{Int64}}` valence-core
- `qc::Dict{Int64,Vector{Int64}}` qspace-core
- `vv::Dict{Int64,Vector{Int64}}` valence-valence
- `qv::Dict{Int64,Vector{Int64}}` qspace-valence
- `qq::Dict{Int64,Vector{Int64}}` qspace-qspace
"""
mutable struct space_channel
    pp::Dict{Int64,Vector{Int64}}
    ph::Dict{Int64,Vector{Int64}}
    hh::Dict{Int64,Vector{Int64}}
    cc::Dict{Int64,Vector{Int64}}
    vc::Dict{Int64,Vector{Int64}}
    qc::Dict{Int64,Vector{Int64}}
    vv::Dict{Int64,Vector{Int64}}
    qv::Dict{Int64,Vector{Int64}}
    qq::Dict{Int64,Vector{Int64}}
end

"""
mutable struct `SingleParticleState`
# Fields
- `n::Int64` principal quantum number of the single particle state(sps)
- `l::Int64` azimuthal quantum number of the sps
- `j::Int64` angular momentum
- `tz::Int64` z-component of isospin (doubled) tz=-1 => proton & tz=1 => neutron
- `occ::Float64` occupation number (can be fractional) of the sps
- `c::Bool` indicating whether the single-particle state belongs to "core" or not 
- `v::Bool` whether belongs to "valence" or not 
- `q::Bool` whether belongs to "q-space" or not 
"""
mutable struct SingleParticleState
    n::Int64
    l::Int64
    j::Int64
    tz::Int64
    occ::Float64
    c::Bool
    v::Bool
    q::Bool
end

"""
struct `ModelSpace`
# Fields
- `p_sps::Vector{SingleParticleState}` proton single particle states (only odd integer)
- `n_sps::Vector{SingleParticleState}` neutron single particle states
- `sps::Vector{SingleParticleState}` single particle states (odd number ones=>proton even=>neutron)
- `occ_p::Matrix{Float64}` matrix representing the occupation number of proton (needed for density matrix)
- `occ_n::Matrix{Float64}` matrix representing the occupation number of neutron
- `holes::Vector{Vector{Int64}}` idx list of holes
- `particles::Vector{Vector{Int64}}` idx list of particles
- `spaces::space_channel` space_channel (mutable struct)
"""
struct ModelSpace
    p_sps::Vector{SingleParticleState}
    n_sps::Vector{SingleParticleState}
    sps::Vector{SingleParticleState}  
    occ_p::Matrix{Float64}
    occ_n::Matrix{Float64}
    holes::Vector{Vector{Int64}}
    particles::Vector{Vector{Int64}}
    spaces::space_channel
end

"""
mutable struct `Operator`
# Fields
- `zerobody::Vector{Float64}` zerobody part of the operator
- `onebody::Vector{Matrix{Float64}}` one-body matrix elements ([1]=>proton, [2]=>neutron)
- `twobody::Vector{Matrix{Float64}}` two-body matrix elements, having array structure [ch]
- `hermite::Bool` whether it is hermitian operator or not
- `antihermite::Bool` antihermitian or not
"""
mutable struct Operator
    zerobody::Vector{Float64}
    onebody::Vector{Matrix{Float64}}
    twobody::Vector{Matrix{Float64}}
    hermite::Bool
    antihermite::Bool
end

"""
struct `HamiltonianNormalOrdered`
# Fields
- `H::Operator` Hamiltonian operator
- `E0::Float64` NO0B of H
- `EMP2::Float64` PT2 correction to E0
- `EMP3::Float64` PT3 correction to E0
- `Cp::Matrix{Float64}` eigenvectors of hp, used unitary trans. HO=>HF basis (proton part)
- `Cn::Matrix{Float64}` eigenvectors of hn, unitary trans. HO=>HF basis (neutron part)
- `e1b_p::Vector{Float64}` eigenvalues of hp
- `e1b_n::Vector{Float64}` eigenvalues of hn
- `modelspace::ModelSpace`
"""
struct HamiltonianNormalOrdered
    H::Operator
    E0::Float64
    EMP2::Float64
    EMP3::Float64
    Cp::Matrix{Float64}
    Cn::Matrix{Float64}
    e1b_p::Vector{Float64}
    e1b_n::Vector{Float64}
    modelspace::ModelSpace
end

"""
mutable struct `IMSRGobject`
# Fields
- `H0::Operator` Hamiltonian for starting point of BCH product
- `H::Operator` Hamiltonian ``H(s)``
- `s::Vector{Float}` current ``s`` and ``ds``
- `smax::Float` maximum ``s``
- `dsmax::Float` maximum ``ds``
- `maxnormOmega::Float` maximum ||Omega||
- `eta::Operator` generator of IMSRG flow (antihermite Operator)
- `Omega::Operator` generator of IMSRG flow (antihermite Operator) 
- `eta_criterion::Float` ||eta|| to check convergence
- `denominatorDelta::Float64` parameter for multi-major shell decoupling
- `n_written_omega::Int` # of written Omega by splitting to solve IMSRGflow
- `Ncomm::Vector{Int}` # of commutator evaluated during IMSRG flow
"""
mutable struct IMSRGobject
    H0::Operator
    H::Operator
    s::Vector{Float64} # [s,ds]
    smax::Float64
    dsmax::Float64
    maxnormOmega::Float64
    eta::Operator
    Omega::Operator
    eta_criterion::Float64
    denominatorDelta::Float64
    n_written_omega::Vector{Int64}
    Ncomm::Vector{Int64}
end

"""
struct `dWS2n`, Wigner symbols used in PreCalcHOB
# Fields
- `dtri::Dict{Vector{Int64},Float64}` dict for trinomial 
- `dcgm0::Dict{Vector{Int64},Float64}` dict for special CG coefficients
- `keycg::Vector{Vector{Int64}}` array of key for cg
"""
struct dWS2n
    dtri::Dict{Vector{Int64},Float64}
    dcgm0::Dict{Vector{Int64},Float64}
    keycg::Vector{Vector{Int64}}
end
