# IMSRG

Files in `src/IMSRG.jl` for IM-SRG calculations:
- [imsrg_util.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/IMSRG.jl/imsrg_util.jl): contains main and util functions 
- [commutator.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/IMSRG.jl/commutator.jl): functions to calculate commutators and BCH transform to carry out IMSRG flow with Magnus expansion
- [valencespace.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/IMSRG.jl/valencespace.jl): functions for Valence-space IM-SRG (VS-IMSRG) calculations to derive shell-model effective interactions/operators


```@autodocs
Modules = [NuclearToolkit]
Pages = ["IMSRG.jl/commutator.jl",
         "IMSRG.jl/imsrg_util.jl",
         "IMSRG.jl/valencespace.jl"]
``` 
