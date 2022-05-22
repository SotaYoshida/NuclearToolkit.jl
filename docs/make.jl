using NuclearToolkit
using Documenter

DocMeta.setdocmeta!(NuclearToolkit, :DocTestSetup, :(using NuclearToolkit); recursive=true)

makedocs(;
    modules=[NuclearToolkit],
    authors="SotaYoshida <s.yoshida@nt.phys.s.u-tokyo.ac.jp> and contributors",
    repo="https://github.com/SotaYoshida/NuclearToolkit.jl/blob/{commit}{path}#{line}",
    sitename="NuclearToolkit.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://SotaYoshida.github.io/NuclearToolkit.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SotaYoshida/NuclearToolkit.jl",
    devbranch="main",
)
