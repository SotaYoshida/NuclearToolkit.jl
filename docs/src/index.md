# NuclearToolkit

Julia Toolkit for nuclear structure calculations

## Installation and example

First, prepare Julia environment v >= 1.7.0.  

Second, add the package in Pkg mode
```julia
julia>]add NuclearToolkit
```

Too see how to run the code, it is strongly recommended to clone the repository.

sample script `sample_script.jl` in the repository performs
- calculating NN potential by Chiral EFT
- HFMBPT(3) and IMSRG/VS-IMSRG(2) calculation
- shell-model calculations with the effective interaction derived by VS-IMSRG

in sequence. One can try as follows
```bash
julia -t 8 sample_script.jl
```

You can specify target nuclei, emax, hw, etc., by editting the script and `optional_parameters.jl`, if necessary.

## Package features and building blocks

NuclearToolkit.jl provides a self-contained set of nuclear structure calculation codes covering from nuclear force to many-body methods (HF/HFMBPT, IM-SRG/VS-IMSRG, shell model, etc.).

One of the main motivations for the author to develop NuclearToolkit.jl is, of course, for their own research purposes, and another one is for educational purpose.
A single nuclear structure calculation code often contains tens of thousands to nearly 100,000 lines.
In the author's personal opinion, it would be too difficult for students (especially in undergraduate or two-year master course) to understand the technical details of the latest nuclear many-body methods while reading the enormous length of existing codes in the community.

The author thought the Julia language can be a game changer to this situation with its high readbility, portabillity, and performance. Since all the source code in NuclearToolkit.jl is a pure Julia implementation, there is no more need to prepare different Makefiles for different environments, worry about library dependencies, homemade Python script to run the Fortran/C++ codes. The code can be easily executed on either a laptop or a supercomputer.
While NuclearToolkit covers a wide range of methods, the overall code length is still in a few tens of thousands, including "docstring" to generate the document.

- [ChiEFTint](ChiEFTint.md): NN interaction from Chiral EFT ~ 6,000 lines.
  - EM (Entem & Machleidt) N3L0
  - EMN (Entem, Machleidt, Nosyc) N4LO
  - Density-Dependent NN from 3NF (The author prefers to call it "2n3n" to distinguish with genuine 3NF)
  - valence chiral EFT potential upto LO
- [Hartreefock](HartreeFock.md): Hartree-Fock (HF) and HF Many-Body Perturbation Theory (HFMBPT)  ~ 3,000 lines.
  - Energy (up to 3rd order)
  - Scaler operator (up to 2nd order)
- [IMSRG](IMSRG.md): In-medium Similarity Renormalization Group (IMSRG)  ~ 2,000 lines.
  - IMSRG(2) calc. for ground state energy
  - consistent IMSRG(2) flow of (scaler) operator
  - Valence-space IMSRG (VS-IMSRG)
    - derive effective interaction for shell-model calculations
    - consistent VSIMSRG flow to get effective operators 
- [ShellModel.jl](ShellModel.md) ~ 5,000 lines.
  This was originally developed as [an independent package] (https://github.com/SotaYoshida/ShellModel.jl).
  - shell model calculations
  - construct approximate wavefunctions with eigenvector continuation 

## Optional parameters
For some parameters, the default values are used unless the user specifies those in the file, `optional_parameters.jl`.
See the [Optional parameters](parameters.md) page for more details.

## Issues/Pull requests

NuclearToolkit.jl is designed to be an open-source software and to guarantee reproducibility and transparancy of the future works.
Making issues and pull requests are fully welcome.
