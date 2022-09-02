# snt format

In NuclearToolkit.jl, input/output interactions are in the so-called snt format (or its binary form).
The snt format is one used in widely-used shell-model code [KSHELL](https://sites.google.com/alumni.tsukuba.ac.jp/kshell-nuclear/).
If you specify the option `tbme_fmt="snt.bin"` for Chiral EFT interactions, the output becomes binary form, where the integers are in `Int64` and floating point numbers are in `Float64`.

Here is an example of w.snt. The comments are not written in the outputs of NuclearToolkit.
```
! "!" denote comment line
!
! model space
! proton-orbit, neutron-orbit, proton core, neutron core
3 3 8 8
! orbit n l j tz
1 0 2 3 -1 ! 1 = p 0d_3/2
2 0 2 5 -1 ! 2 = p 0d_5/2
3 1 0 1 -1 ! 3 = p 1s_1/2
4 0 2 3 1 ! 4 = n 0d_3/2
5 0 2 5 1 ! 5 = n 0d_5/2
6 1 0 1 1 ! 6 = n 1s_1/2
! one-body interaction
! number of lines, method1
6 0
! i j <i|V|j>
1 1 1.64658
2 2 -3.94780
3 3 -3.16354
4 4 1.64658
5 5 -3.94780
6 6 -3.16354
! two-body interaction (TBME)
! # of lines, method2 A mass dependence factor
158 1 18 -0.30000
! i j k l J <i,j| V | k,l>_J
1 1 1 1 0 -2.18450
1 1 1 1 2 -0.06650
```
