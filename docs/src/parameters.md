# Parameters in NuclearToolkit

Here is the summary of optional parameters, which can be specified by users through `optional_parameters.jl` in the working directory.

* For ChiEFTint
    - `n_mesh::Int64`, number of momentum mesh 
    - `pmax_fm::Float64`, maximum pmax in fm``^{-1}``
    - `emax::Int64`, emax quanta
    - `Nnmax::Int64`, Nnmax quanta
    - `chi_order::Int64`, order of chiral EFT potenteial (0:LO 1:NLO 2:NNLO 3:N3LO 4:N4LO)
    - `calc_NN::Bool`, calculate NN part
    - `calc_3N::Bool`, calculate eff2n (2n3n) part
    - `hw::Float64` hbaromega
    - `srg::Bool`, SRG transformation (NN-only)
    - `srg_lambda::Float64`, resolution scale for SRG in ``fm^{-1}``
    - `tbme_fmt::String` "snt" or "snt.bin" is supported
    - `fn_tbme::String` file name of output interaction
    - `pottype::String` potential type (em500n3lo,emn500n3lo,emn500n4lo) 
    - `kF::Float64` Fermi momentum for 2n3n
* For IMSRG
    - `smax::Float64` max flow parameter `s`
    - `dsmax::Float64` maximum step size for IMSRG flow
    - `maxnormOmega::Float64` tol for norm of IMSRG generator
    - `denominatorDelta::Float64` denominator delta for multi-shell interaction