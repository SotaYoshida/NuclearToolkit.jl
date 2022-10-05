# Contributing to NuclearToolkit.jl

Thank you for considering contributing to this package.

Feedbacks and contributions to NuclearToolkit.jl are very welcome.
These can be:
- bug report
- submitting a new function or a patch to the bug
- documentation issue
- feature request
- etc.

For these contributions, it would be nice to let you know a basic policy (the workflow of code development and LICENSE, shown below) in this package.
Note that the package has been developed by a single author (@SotaYoshida) so far, and thereby the followings are just the author's policy. Comments on the development policy itself are also welcome.


## Workfklow of code development

We use the GitHub to host the package, to track issues/pull requests.

The document is built using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/),
which is a package for building documentation from docstrings and markdown files.
It is automized to build and deploy the document by GitHub Actions.

A typical workfklow is the following:

1. clone or folk the repository and modify the code and/or document locally

2. propose modifications through a pull request (PR) to **'dev' branch**

3. if the PR passes the automated tests, it will be merged 

4. At some point, we tag and release, which is done automatically by instructing the JuliaRegistrator bot to do so.


### Automated tests on GitHub Actions

We use GitHub Actions to run the test codes and to build/deploy the document.
When some changes are submitted through a pull request, the test codes are run to check that the changes are not destructive.

The test jobs are specified in yml files like `.github/workflows/CI.yml` and one can find the test code in `test/` of the repository.

They checks
- Generating a nucleon-nucleon potential, i.e., input for nuclear many-body methods
- HFMBPT results with the NN-potential
- IMSRG/VSIMSRG results the NN-potential
- shell-model results with the derived effective interaction by VS-IMSRG
- MPI and sampling stuffs are runnable

If you submit a major change to the code and it brings the changes in the results by many-body methods,
it is recommended to contact the main author(s). Otherwise, we appreciate it if you could propose updated test codes.


## Reporting bugs by opening a new issue

Thank you so much for considering to report bugs!
When you report a bug in the code, please open an new issue from [Here](https://github.com/SotaYoshida/NuclearToolkit.jl/issues).
Please make sure to include any information necessary for us to reproduce the errors. Thanks!

## Propose modifications through Pull Requests (PRs)

You can propose modifications you made on the package through pull requests.

* As stated above, please consider to make test codes matching to your modifications.
 
* Please make sure to submit your PR to `develop` branch. The 'main' branch will be protected by 'github branch protection'.

As this package is currently being developed by a single author, branching rules such as git-flow and GitHub Flow have not been adopted.
When we got contributors, the branch-rule will be set upon discussions.


### Style Guide for Julia

If you are new to Julia language, Julia's [Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/) can be usefull to write your own codes. If the current code does not follow the style guide, we would be grateful if you could let us know.

It is recommended (not always necessarily) to add the so-called 'docstring' to the functions to explain what the function does and to generate document with Documenter.jl. An example is:
```julia
"""
  Eval_EigenStates!(H::Hamiltonian,ψ::WaveFunction,n::Int64)

Function to calculate the lowest ``n`` wavefunctions and overwrite results to `ψ`.

# Arguments
- `H::Hamiltonian` the Hamiltonian of the system, see struct `Hamiltonian` for more details.
- `ψ::WaveFunction` the wavefunction, see struct `WaveFunction` for more details.
- `n::Int64` the number of states of interest.
"""
function Show_EigenStates!(H::Hamiltonian,ψ::WaveFunction,n::Int64)
    # write your nice codes here
    
end
```

When you want to add a new API to the package, which can be called in Julia REPL, you need to export it `src/NuclearToolkit.jl`.


## LICENSE

Any contribution from you will be under the MIT License, as well as the package itself.
Feel free to contact to @SotaYoshida if that's a concern.


