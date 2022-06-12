# NuclearToolkit

Julia Toolkit for nuclear structure calculations

## Installation and example

First, prepare Julia environment v >= 1.7.0.  

Second, add the package in Pkg mode
```julia
julia>]add NuclearToolkit
``` 
Note: The above is currently not working, since this package has not yet registered as an official package.  
The adhoc prescription for now is to execute below
```bash
julia src/package_install.jl
```

A sample script provided to perform
- calculating NN potential by Chiral EFT
- HFMBPT(3) and IMSRG/VS-IMSRG(2) calculation
- shell-model calculations with the effective interaction derived by VS-IMSRG

in sequence. One can try as follows
```bash
julia -t 8 sample_script.jl
```

If you want to specify target nuclei, modelspace, hw, etc., please edit `optional_parameters.jl`.

## Package features and building blocks

NuclearToolkit.jl provides a self-contained set of nuclear structure calculation codes covering from nuclear force to many-body methods (HF/HFMBPT, IM-SRG/VS-IMSRG, shell model, etc.).

One of the main motivations for the author to develop NuclearToolkit.jl is, of course, for their own research purposes, and another one is for educational purpose.
A single nuclear structure calculation code often contains tens of thousands to nearly 100,000 lines.
In the author's personal opinion, it would be too difficult, students (especially in undergraduate or two-year master course), to understand the technical details of the latest nuclear many-body methods while reading the enormous length of existing codes in the community.

The author thought the Julia language can be a game changer to this situation with its high readbility, portabillity, and performance. Since all the source code in NuclearToolkit.jl is a pure Julia implementation, there is no more need to prepare different Makefiles for different environments or worry about library dependencies. One can work equally well on a laptop or a supercomputer.
While NuclearToolkit covers a wide range of methods, the overall code length is still in a few tens of thousands, including "docstring" to generate the document.

- [ChiEFTint](build/ChiEFTint.html): NN interaction from Chiral EFT ~ 5,000 lines.
  - Entem & Machleidt N3L0
  - Density-Dependent NN from 3NF
  - "Valence" NO
- [Hartreefock](../build/HartreeFock.html): Hartree-Fock (HF) and HF Many-Body Perturbation Theory (HFMBPT)  ~ 3,000 lines.
  - Energy (up to 3rd order)
  - Scaler operator (up to 2nd order)
- [IMSRG](build/IMSRG.html): In-medium Similarity Renormalization Group (IMSRG)  ~ 2,000 lines.
  - IMSRG(2) calc. for ground state energy
  - consistent IMSRG(2) flow of (scaler) operator
  - Valence-space IMSRG (VS-IMSRG)
    - derive effective interaction for shell-model calculations
    - consistent VSIMSRG flow to get effective operators 
- [ShellModel.jl](build/ShellModel.html) ~ 5,000 lines.
  This was originally developed as [an independent package] (https://github.com/SotaYoshida/ShellModel.jl).
  - shell model calculations
  - construct approximate wavefunctions with eigenvector continuation 

## Optional parameters
For some parameters, the default values are used unless the user specifies some options in the file `optional_parameters.jl`.
See the [Optional parameters](build/parameters.html) page for more details.

## Issues/Pull requests

NuclearToolkit.jl is designed to be an open-source software
and to guarantee reproducibility and transparancy of future works.
Making issues and pull requests are fully welcome.
