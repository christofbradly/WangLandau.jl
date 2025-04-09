using WangLandau
using Documenter

DocMeta.setdocmeta!(WangLandau, :DocTestSetup, :(using WangLandau); recursive=true)

makedocs(;
    modules=[WangLandau],
    authors="Chris Bradly",
    sitename="WangLandau.jl",
    format=Documenter.HTML(;
        canonical="https://christofbradly.github.io/WangLandau.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Guide" => "guide.md",
        "Example: Ising model" => "ising.md",
        "Advanced usage" => "advanced.md",
        "Index" => "index.md"
    ],
)

deploydocs(;
    repo="github.com/christofbradly/WangLandau.jl",
    devbranch="main",
)
