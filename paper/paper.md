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
date: 10 August 2022
bibliography: paper.bib
---

# Summary
One of the ultimate goals in nuclear physics is to understand and predict various properties of nuclei from a fundamental interaction among nucleons, the nuclear force.
Owing to recent developments in the community, especially in describing the nuclear force and in nuclear many-body methods, it is becoming possible to make quantitative discussions and predictions on various properties with first-principles calculations.
On the other hand, it is a common situation that different degrees of freedom can be important (e.g., nucleons, alpha-cluster) to describe one nucleus or one state of a target nucleus, so nuclear models are becoming more and more diverse.
It can therefore be difficult to familiarize oneself with the technical details of all those methods. ``NuclearToolkit.jl`` was designed to circumvent the situation and to be helpful for both students and researchers to tackle nuclear many-body problems.

# Statement of need

``NuclearToolkit.jl`` provides self-contained codes for nuclear physics covering from nuclear forces to various nuclear many-body methods.
Users can generate nucleon-nucleon (NN) potentials based on chiral effective field theory [@EMrev1;@EMrev2;@LENPICrev1;@LENPICrev2] and use them in many-body methods such as Hartree-Fock many-body perturbation theory (which is well known as Møller–Plesset method in the community of chemistry) [@ShavittBartlett2009], in-medium similarity renormalization group (IM-SRG) [@StrobergRev19], and valence shell model (configuration interaction method) [@BrownRev;@CaurierRev;@OtsukaRev].

In nuclear physics community, many public codes are already available.
Although it may not be possible to list them all, representative examples are ANTOINE [@ANTOINE], NuShellX [@NuShellX], BIGSTICK [@BIGSTICK1;@BIGSTICK2], and KSHELL [@KSHELL1;@KSHELL2] for valence shell-model and imsrg [@imsrgcode] for IM-SRG.
People are using those various codes, which are typically written in Fortran, C++, etc., and call them in their homemade shell or Python scripts.
`NuclearToolkit.jl`, which is not a wrapper of those existing codes, provides
a new interface that combines these various methods into one and works on a variety of environments, including Linux, Mac, and Windows.
This is achieved thanks to the high readability and portability of the Julia programming language[@Bezanson2012].

This code has been used in previous works [@PTEPSY;@SY2n3n].

# Acknowledgements

The author thanks Noritaka Shimizu, Takayuki Miyagi, and Tokuro Fukui 
for discussions at the early stage of the development of NuclearToolkit.jl,
and Michio Kohno for discussions on the density-dependent NN force.
This work was supported by JSPS KAKENHI (Grants No. 22K14030) and partially by the computational resource of Fujitsu PRIMERGY CX400M1/CX2550M5 (Oakbridge-CX) at the Information Technology Center, The University of Tokyo.


# References