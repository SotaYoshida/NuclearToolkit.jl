# NuclearToolkit

Julia Toolkit for nuclear structure calculations

## Installation and example

First, prepare Julia environment v >= 1.7.0.  

Second, add the package in Pkg mode
```julia
julia>]add NuclearToolkit
```

You can try the package in the following ways:

* clone the repository and run `example/sample_script.jl` in the repository like

    `$ julia -t 8 example/sample_script.jl`

    This performs:
      - calculation of NN potential from Chiral EFT
      - HFMBPT(3) and IMSRG/VS-IMSRG(2) calculation with it
      - shell-model calculations with the effective interaction derived by VS-IMSRG
    An expected results using the latest dev branch can be found [here](https://github.com/SotaYoshida/NuclearToolkit.jl/blob/develop/example/log_sample_script.txt).
* Try sample codes in [HowToUse](howtouse) page.



Please make sure to use the latest version of the package. Update can be done with 
```julia
julia>]up NuclearToolkit
```
In the Julia REPL, you can see the UUIDs and versions of the installed packages
```julia
julia>using Pkg
julia>Pkg.status()
```

## Package features 

NuclearToolkit.jl provides a self-contained set of nuclear structure calculation codes covering from nuclear force to many-body methods (HF/HFMBPT, IM-SRG/VS-IMSRG, shell model, etc.).

One of the main motivations for the author to develop NuclearToolkit.jl is, of course, for their own research purposes, and another one is for educational purpose.
A single nuclear structure calculation code often contains tens of thousands to nearly 100,000 lines.
In the author's personal opinion, it would be too difficult for students (especially in undergraduate or two-year master course) to understand the technical details of the latest nuclear many-body methods while reading the existing code with enormous length.

The author thought the Julia language can be a game changer to this situation with its high readbility, portabillity, and performance. Since all the source code in NuclearToolkit.jl is a pure Julia implementation, there is no more need to prepare different Makefiles for different environments, worry about library dependencies, homemade Python scripts to run the Fortran/C++ codes. The code can be easily executed on different environments, including Linux, Mac, and Windows.
While NuclearToolkit covers a wide range of methods, the overall code length is still in a few tens of thousands, including "docstring" to generate the document.

## Building blocks

- [ChiEFTint](ChiEFTint): NN interaction from Chiral EFT ~ 6,500 lines.
  - EM (Entem & Machleidt) N3L0
  - EMN (Entem, Machleidt, Nosyk) N4LO
  - Density-Dependent NN from 3NF (The author prefers to call it "2n3n" to distinguish with genuine 3NF)
  - valence chiral EFT potential upto LO
- [Hartreefock](HartreeFock): Hartree-Fock (HF) and HF Many-Body Perturbation Theory (HFMBPT)  ~ 3,000 lines.
  - Energy (up to 3rd order)
  - Scaler operator (up to 2nd order)
- [IMSRG](IMSRG): In-medium Similarity Renormalization Group (IMSRG)  ~ 2,000 lines.
  - IMSRG(2) calc. for ground state energy
  - consistent IMSRG(2) flow of (scaler) operator
  - Valence-space IMSRG (VS-IMSRG)
    - derive effective interaction for shell-model calculations
    - consistent VSIMSRG flow to get effective operators 
- [ShellModel.jl](ShellModel): valence shell model code ~ 5,000 lines.  

  This was originally developed as [an independent package] (https://github.com/SotaYoshida/ShellModel.jl).
  - shell model calculations
  - construct approximate wavefunctions with eigenvector continuation 

## Optional parameters
For some parameters, the default values are used unless the user specifies.
See the [Optional parameters](parameters) page for more details.

## Issues/Pull requests

NuclearToolkit.jl is designed to be an open-source software and to guarantee reproducibility and transparancy of the future works.
Making issues and pull requests are fully welcome.
If you would like to contribute to development of the package, see the [Contributing](contributing) page and please do not hesitate to contact me (S.Yoshida).

