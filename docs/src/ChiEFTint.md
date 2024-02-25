# ChiEFTint

Functions needed to generate nucleon-nucleon (NN) potential from Chiral EFT.

The parameters needed for chiEFTint are specified through `optional_parameters.jl` or the optional argument `fn_params=[PATH_TO_FILE]` in main API, `make_chiEFTint()`.

```@autodocs
Modules = [NuclearToolkit]
Pages = ["chiEFTint/main_chiEFTint.jl",
         "chiEFTint/contact.jl",
         "chiEFTint/pionexchange.jl",
         "chiEFTint/angmom_algebra.jl",
         "chiEFTint/misc_plt_io.jl",
         "chiEFTint/threebodyforce.jl",
         "chiEFTint/calibration.jl",
         "chiEFTint/renorm.jl",
         "chiEFTint/valence.jl",         
         "chiEFTint/eff3nf.jl",
	     "chiEFTint/read_me3j.jl",
         "chiEFTint/dict_LECs.jl"
]
``` 

