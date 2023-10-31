# IMSRG

Files in `src/IMSRG.jl` for IM-SRG calculations:
- [imsrg_util.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/IMSRG.jl/imsrg_util.jl): contains main and util functions 
- [commutator.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/IMSRG.jl/commutator.jl): functions to calculate commutators and BCH transform to carry out IMSRG flow with Magnus expansion
- [valencespace.jl](https://github.com/SotaYoshida/NuclearToolkit.jl/tree/main/src/IMSRG.jl/valencespace.jl): functions for Valence-space IM-SRG (VS-IMSRG) calculations to derive shell-model effective interactions/operators



Since we use the so-called Magnus formulation of IMSRG flow, $\Omega$

```math
H(s) = e^{\Omega(s)}H(0)e^{-\Omega(s)} \\\\
e^{\Omega(s+ds)} = e^{\eta(s)ds}e^{\Omega(s)}
```
is needed to explicitly calculate the flow generator ``\eta``.
See e.g. [Annual Review of Nuclear and Particle Science 69, 307-362](https://doi.org/10.1146/annurev-nucl-101917-021120) for more details.

In the package, temporary binary files for ``\Omega`` are genereted in `flowOmega` directory (which will be made if not available).
Although users do not need to read or write these binary files in normal use, the contents of the binary files are described for specific uses and development, e.g. re-start IMSRG flow from the temporary files.

A temporary file is generated like `flowOmega/Omega_1234He4_1.bin`.
This contains
- number of single particle basis (`Int64`), which is assumbed to be common between proton and neutron.
- One-body matrix element for protons and neutrons (`Matrix{Float64}`)
- number of two-body channels (`Int64`)
- dimensions of each channel (`Vector{Int64}`)
- Two-body matrix element for protons and neutrons (`Vector{Matrix{Float64}}`)

Note that the all matrix elements are written in row-major order (I know Julia is column-major language).

---

```@autodocs
Modules = [NuclearToolkit]
Pages = ["IMSRG.jl/commutator.jl",
         "IMSRG.jl/imsrg_util.jl",
         "IMSRG.jl/valencespace.jl"]
``` 
