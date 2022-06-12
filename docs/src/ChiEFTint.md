# ChiEFTint

Functions needed to generate nucleon-nucleon (NN) potential from Chiral EFT.

The parameters needed for chiEFTint are specified in `optional_parameters.jl`.
If one specifies the potential type in `optional_paramters.jl`, the corresponding `LECs_XX.jl` will be read (if possible).
Otherwise, one runs with the so-called Entem Machleidt N3LO potential.

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

