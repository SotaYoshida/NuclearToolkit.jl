# NuclearToolkit.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://SotaYoshida.github.io/NuclearToolkit.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://SotaYoshida.github.io/NuclearToolkit.jl/dev)
[![Build Status](https://github.com/SotaYoshida/NuclearToolkit.jl/actions/workflows/CI.yml/badge.svg?branch=dev)](https://github.com/SotaYoshida/NuclearToolkit.jl/actions/workflows/CI.yml?branch=dev)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04694/status.svg)](https://doi.org/10.21105/joss.04694)
<img src="https://github.com/SotaYoshida/NuclearToolkit.jl/blob/main/docs/src/assets/logo_full.png" width=60%>


Julia package for nuclear structure calculations covering:
- generating Chiral EFT interactions
- many-body calculations (HFMBPT, IMSRG/VS-IMSRG, valence shell-model, etc.)


Note: Of course, 'for structural calculations' simply means that the author (SY) is not familiar with reaction theories.
Contributions and suggestions from reaction theory and experimental researchers are very welcome. Thanks.

## Installation

Assuming that you have already installed Julia (v>=1.7.0),
```jldoctest
julia>using Pkg; Pkg.add("NuclearToolkit")
```
 
Pkg (Julia's builtin package manager) intalls packages in ```$JULIA_DEPOT_PATH```, which is by default ```~/.julia```.  
When working on a working node (w/o permissions to access ```~/```), overwrite the ```JULIA_DEPOT_PATH``` by ```export JULIA_DEPOT_PATH="PATH_TO_JULIA_DEPOT"```.

## How to START

Execute `example/sample_script.jl` like
```
julia -t 10 example/sample_script.jl
```

This sample script performs:
 - generating NN potential with hw=20, emax=4
 - HFMBPT(3) and IMSRG/VS-IMSRG calculation using the NN potential 
 - shell-model calculation using the effective interaction derived by VS-IMSRG

## How to cite

When you use `NuclearToolkit.jl` in your work, please cite the JOSS paper:https://doi.org/10.21105/joss.04694
```
@article{Yoshida2022,
doi = {10.21105/joss.04694}, 
url = {https://doi.org/10.21105/joss.04694}, 
year = {2022}, 
publisher = {The Open Journal}, 
volume = {7}, number = {79}, pages = {4694}, 
author = {Sota Yoshida}, 
title = {NuclearToolkit.jl: A Julia package for nuclear structure calculations}, 
journal = {Journal of Open Source Software} }
```