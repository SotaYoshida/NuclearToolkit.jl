---
title: 'NuclearToolkit.jl: A Julia package for nuclear structure calculations'
tags:
  - Julia
  - nuclear physics
  - nuclear force
  - nuclear structure
  - ab initio calculation
authors:
  - name: Sota Yoshida
    orcid: 0000-0002-1342-1846
    affiliation: 1
affiliations:
 - name: Institute for Promotion of Higher Academic Education, Utsunomiya University
   index: 1
date: 16 May 2022
bibliography: paper.bib
---


# Summary
One of the ultimate goals in nuclear physics is to understand and predict various properties of nuclei from a fundamental interaction among nucleons, nuclear force.
Owing to recent developments in the community, especialy in describing the nuclear force and in nuclear many-body methods, it is becoming possible to make quantitative discussions and predictions on various properties with first-principles calculations.
On the other hand, it is a common situation that different degrees of freedom can be important (e.g., nucleons, alpha-cluster) to describe one nucleus or one state of a target nucleus, so nuclear models are becoming more and more diverse.
It can therefore be difficult to familiarize oneself with the technical details of all those methods. ``NuclearToolkit.jl`` was designed to circumvent the situation and to be helpful for both students and researchers to tackle nuclear many-body problems.

# Statement of need

``NuclearToolkit.jl`` provides self-contained codes covering from nuclear forces to various nuclear many-body methods.
Users can generate nucleon-nucleon (NN) potentials based on chiral effective field theory [@RevChiEFT_Idaho,@RevChiEFT_LENPIC] and use them in many-body methods such as Hartree-Fock many-body perturbation theory (which is well known as Møller–Plesset method in the community of chemistry) [@ShavittBartlett2009], in-medium similarity renormalization group (IM-SRG) [@StrobergRev19], and valence shell model (configuration interaction method) [@BrownRev,@CourierRev,@OtsukaRev].

In nuclear physics community, many public codes are already available.
Although it may not be possible to list them all, ANTOINE [@ANTOINE], NuShellX [@NuShellX], BIGSTICK [@BIGSTICK1,@BIGSTICK2], and KSHELL [@KSHELL1,@KSHELL2] for valence shell-model and imsrg [@imsrgcode] for IM-SRG.
People are using those various codes, which are typically written in FORTRAN77, Fortran90, C++, etc., and call them in their homemade codes.
`NuclearToolkit.jl`, which is not a wrapper of those existing codes, provides
a new interface that combines these various methods into one and works on both laptops and supercomputers in almost the same way. That is achieved by high readability and portability of Julia language[@Bezanson2012].

A part of the code, which was originally developed as an independent package ``ShellModel.jl`` [@ShellModel.jl], have already been used in a published work [@PTEP_SY], and the IM-SRG method implemented in the package is able to derive shell-model effective interactions for the model space of interest and is becoming a kind of gold standard for studying open-shell nuclei from first principles.
Therefore, the author expects that many future studies will emerge to tackle the cutting edge of nuclear many-body problems with this package.


# Acknowledgements

This work was supported by JSPS KAKENHI Grant Number 22K14030.
The author acknowledges Noritaka Shimizu and Takayuki Miyagi for discussions at early stage of the code development.

# References