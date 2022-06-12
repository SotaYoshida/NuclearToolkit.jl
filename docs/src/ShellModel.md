# ShellModel

Functions for shell-model calculations
- [shellmodel_main.jl](../../src/ShellModel/shellmodel_main.jl): contains main and util functions 
- [lanczos_methods.jl](../../src/ShellModel/lanczos_methods.jl): Lanczos methods
- [transit.jl](../../src/ShellModel/transit.jl): EM transitions
- [input_int_snt.jl](../../src/ShellModel/input_int_snt.jl): I/O stuffs
- [eigenvector_continuation.jl](../../src/ShellModel/eigenvector_continuation.jl): eigenvector continuation to sample and construct shell-model wavefunctions


```@autodocs
Modules = [NuclearToolkit]
Pages = ["ShellModel/shellmodel_main.jl",
         "ShellModel/lanczos_methods.jl",
         "ShellModel/transit.jl",
         "ShellModel/input_int_snt.jl",
         "ShellModel/eigenvector_continuation.jl"]

``` 
