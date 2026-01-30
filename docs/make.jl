using ARDmodels
using Documenter

DocMeta.setdocmeta!(ARDmodels, :DocTestSetup, :(using ARDmodels); recursive=true)

makedocs(;
    modules=[ARDmodels],
    authors="Laurent Doyen",
    sitename="ARDmodels.jl",
    format=Documenter.HTML(;
        canonical="https://ldoyen.github.io/ARDmodels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ldoyen/ARDmodels.jl",
    devbranch="main",
)
