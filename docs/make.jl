using Algoim
using Documenter

DocMeta.setdocmeta!(Algoim, :DocTestSetup, :(using Algoim); recursive=true)

makedocs(;
    modules=[Algoim],
    authors="Eric Neiva <eric.neiva@college-de-france.fr>",
    repo="https://github.com/ericneiva/Algoim.jl/blob/{commit}{path}#{line}",
    sitename="Algoim.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ericneiva.github.io/Algoim.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ericneiva/Algoim.jl",
    devbranch="main",
)
