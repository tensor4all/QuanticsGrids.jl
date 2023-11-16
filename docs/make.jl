using QuanticsGrids
using Documenter

DocMeta.setdocmeta!(QuanticsGrids, :DocTestSetup, :(using QuanticsGrids); recursive=true)

makedocs(;
    modules=[QuanticsGrids],
    authors="Ritter.Marc <Ritter.Marc@physik.uni-muenchen.de> and contributors",
    repo="https://gitlab.com/tensors4fields/QuanticsGrids.jl/blob/{commit}{path}#{line}",
    sitename="QuanticsGrids.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        repolink="https://gitlab.com/tensors4fields/QuanticsGrids.jl",
        edit_link="main",
        assets=String[]),
    pages=[
        "Home" => "index.md",
        "API Reference" => "apireference.md",
    ])
