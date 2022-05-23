# NuclearToolkit.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://SotaYoshida.github.io/NuclearToolkit.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://SotaYoshida.github.io/NuclearToolkit.jl/dev)
[![Build Status](https://github.com/SotaYoshida/NuclearToolkit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SotaYoshida/NuclearToolkit.jl/actions/workflows/CI.yml?query=branch%3Amain)

<img src="https://github.com/SotaYoshida/NuclearToolkit.jl/blob/main/docs/src/assets/logo_full.png" width=60%>


Julia Toolkit for nuclear structure calculations covering,
- generating Chiral EFT interactions
- many-body calculations (HFMBPT, IMSRG/VS-IMSRG, valence shell-model, etc.)

## Installation

Assuming that you have already installed Julia, execute below (these are identical) 
```jldoctest
julia>import Pkg; Pkg.add("NuclearToolkit")
```
or 
```jldoctest
julia>]add NuclearToolkit.jl
```
 
Pkg (Julia's builtin package manager) intalls packages in ```$JULIA_DEPOT_PATH```, which is by default ```~/.julia```.  
When working on a working node (w/o permissions to access ```~/```), overwrite the ```JULIA_DEPOT_PATH``` by ```export JULIA_DEPOT_PATH="PATH_TO_JULIA_DEPOT"```.

## How to START

Execute `sample_script.jl` like
```
julia -t 6 sample_script.jl
```

This sample script do 
 - generate NN potential with hw=20, emax=4
 - HFMBPT(3) and IMSRG/VS-IMSRG calculation using the NN potential 
 - shell-model calculation using the effective interaction derived by V-SIMSRG
