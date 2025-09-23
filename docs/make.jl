using WangLandau
using Documenter

makedocs(;
    modules=[WangLandau],
    authors="Chris Bradly",
    sitename="WangLandau.jl",
    # format=Documenter.HTML(;
    #     canonical="https://christofbradly.github.io/WangLandau.jl",
    #     edit_link="main",
    #     assets=String[],
    #     size_threshold=500_000, # 500 kB
    #     size_threshold_warn=200_000, # 200 kB
    # ),
    pages=[
        "Home" => "index.md",
        "Example: 2D Ising model" => "ising.md",
        "Advanced usage" => "advanced.md",
        "API" => "api.md"
    ],
    # checkdocs=:exports,
    # doctest=false, # Doctests are done while testing.
)

deploydocs(
    repo="github.com/christofbradly/WangLandau.jl"
)
