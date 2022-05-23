using NuclearToolkit
using Documenter
#include("../src/NuclearToolkit.jl")
#push!(LOAD_PATH,"../src/")

#DocMeta.setdocmeta!(NuclearToolkit, :DocTestSetup, :(using NuclearToolkit); recursive=true)
makedocs(;
        modules=[NuclearToolkit],
        authors="SotaYoshida <syoshida@cc.utsunomiya-u.ac.jp>",
        repo="https://github.com/SotaYoshida/NuclearToolkit.jl/blob/{commit}{path}#{line}",
        sitename="NuclearToolkit.jl",
        format=Documenter.HTML(;
            prettyurls=get(ENV, "CI", "false") == "true",
            canonical="https://SotaYoshida.github.io/NuclearToolkit.jl",
            assets=String[],
        ),
    pages=[
        "Home" => "index.md",
        "ChiEFTint" => "ChiEFTint.md",
        "HartreeFock" => "HartreeFock.md",
        "IMSRG" => "IMSRG.md",
        "ShellModel" => "ShellModel.md"
    ],
)

# deploydocs(;
#     repo="github.com/SotaYoshida/NuclearToolkit.jl",
#     devbranch="main",
# )
