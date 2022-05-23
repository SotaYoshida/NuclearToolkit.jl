# ChiEFTint

Functions needed to generate nucleon-nucleon (NN) potential from Chiral EFT.

The parameters needed for chiEFTint are specified in `optional_parameters.jl`.
If you want to change Low-energy constatants, edit the [LECs.jl](../../src/chiEFTint/LECs.jl) or use dictionary for LECs e.g., in [main_chiEFTint.jl](../../src/chiEFTint/main_chiEFTint.jl).

```@autodocs
Modules = [NuclearToolkit]
Pages = ["chiEFTint/main_chiEFTint.jl",
         "chiEFTint/contact.jl",
         "chiEFTint/pionexchange.jl",
         "chiEFTint/angmom_algebra.jl",
         "chiEFTint/misc_plt_io.jl",
         "chiEFTint/threebodyforce.jl",
         "chiEFTint/bayesopt.jl",
         "chiEFTint/valence.jl",         
         "chiEFTint/eff3nf.jl",
]
``` 

