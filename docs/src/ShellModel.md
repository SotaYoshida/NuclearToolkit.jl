# ShellModel

Functions for shell-model calculations
- [shellmodel_main.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/ShellModel/shellmodel_main.jl): contains main and util functions 
- [lanczos_methods.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/ShellModel/lanczos_methods.jl): Lanczos methods
- [transit.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/ShellModel/transit.jl): EM transitions
- [input_int_snt.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/ShellModel/input_int_snt.jl): I/O stuffs
- [eigenvector_continuation.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/ShellModel/eigenvector_continuation.jl): eigenvector continuation to sample and construct shell-model wavefunctions
- [KSHELL.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/ShellModel/KSHELL.jl): dealing with KSHELL outputs
- [betadecay.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/ShellModel/betadecay.jl): evaluate beta-decay properties using KSHELL log/summary

```@autodocs
Modules = [NuclearToolkit]
Pages = ["ShellModel/shellmodel_main.jl",
         "ShellModel/lanczos_methods.jl",
         "ShellModel/transit.jl",
         "ShellModel/input_int_snt.jl",
         "ShellModel/eigenvector_continuation.jl",
         "ShellModel/KSHELL.jl",
         "ShellModel/betadecay.jl"]
``` 
