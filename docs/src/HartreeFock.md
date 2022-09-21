# Hartreefock

Files in `src/hartreefock.jl` for Hartree-Fock calculations.

- [def_struct.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/hartreefock.jl/def_struct.jl): define struct/mutable struct  
- [hf_mbpt.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/hartreefock.jl/hf_mbpt.jl): calculate HFMBPT energy correction
- [io_input.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/hartreefock.jl/io_input.jl): I/O stuffs and read input (snt file)
- [main.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/hartreefock.jl/main.jl): main functions
- [operator.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/hartreefock.jl/operator.jl): (scaler) operators and normal ordering

```@autodocs
Modules = [NuclearToolkit]
Pages = ["hartreefock.jl/def_struct.jl",
        "hartreefock.jl/hf_mbpt.jl",
        "hartreefock.jl/io_input.jl",
        "hartreefock.jl/main.jl",
        "hartreefock.jl/operator.jl"]
``` 
