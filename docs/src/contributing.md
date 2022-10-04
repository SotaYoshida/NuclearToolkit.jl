# Contributing to NuclearToolkit.jl

Feedbacks and contributions to NuclearToolkit.jl are very welcome.
These can be:
- bug report
- submitting a new function or a patch to the bug
- documentation issue
- feature request
- etc.

For these contributions, it would be nice to let you know a basic policy (the workflow of code development, and LICENSE, shown below) in this package.
It should be noted that the package has been developed by a single author (@SotaYoshida) so far, and thereby the followings are just the author's policy. Comments on development policy are also welcome.


## Workfklow of code development

We use the GitHub to host the package, to track issues/pull requests.

We use the Julia package, Documenter.jl, to build the document.
About 
- release
- 

workflow (which branches to submit PR to, etc)

GitHub Actions (.github/workflow/) run the test code

.github/workflows/CI.yml

When you want to add a new API to the package, which can be called in Julia REPL, you need to export it `src/NuclearToolkit.jl`.


### Automated tests on GitHub Actions

We use GitHub Actions to run the test codes and to build/deploy the document.
When some changes are submitted through a pull request, the test codes are run to check that the changes are not destructive.

The test codes and yml files specifying them are in `test` in the repository. They checks
- Generating a nucleon-nucleon potential, i.e., input for nuclear many-body methods
- HFMBPT results with it
- IMSRG/VSIMSRG results with it
- shell-model results with the derived effective interaction
- MPI and sampling stuffs are runnable

If you submit a major change to the code and it brings the changes in the numbers of the result by many-body methods,
it is recommended to contact the main author(s). Otherwise, we appreciate it if you could propose updated test codes.


## Report a bug by opening a new issue

Thank you so much for considering to report bugs!
When you report a bug in the code, please open an issue and make sure to include any information necessary for us to reproduce the bug.


## Propose your modification through Pull Requests (PRs)

You can propose modifications on the package through pull requests.

### Style Guide for Julia

If you are new to Julia language, Julia's [Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/) can be usefull to write your own codes. If the current code does not follow the style guide, we would be grateful if you could let us know.

It is recommended (not always necessarily) to add the so-called 'docstring' to the functions to explain what the function does and to generate document automatically with Documenter.jl. An example is:
```julia
"""
  Eval_EigenStates!(H::Hamiltonian,ψ::WaveFunction,n::Int64)

Function to calculate the lowest ``n`` wavefunctions and overwrite details to `ψ`.

# Arguments
- `H::Hamiltonian` the Hamiltonian of the system, see struct `Hamiltonian` for more details.
- `ψ::WaveFunction` the wavefunction, see struct `WaveFunction` for more details.
- `n::Int64` the number of states of interest.
"""
function Show_EigenStates!(H::Hamiltonian,ψ::WaveFunction,n::Int64)
    # write your nice codes here
    
end
```


## LICENSE

Any contribution from you will be the MIT License, same as the package.
Feel free to contact to @SotaYoshida if that's a concern.


