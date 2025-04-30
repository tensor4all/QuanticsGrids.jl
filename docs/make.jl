using QuanticsGrids
using Documenter

DocMeta.setdocmeta!(QuanticsGrids, :DocTestSetup, :(using QuanticsGrids); recursive = true)

makedocs(;
    modules = [QuanticsGrids],
    authors = "Ritter.Marc <Ritter.Marc@physik.uni-muenchen.de> and contributors",
    sitename = "QuanticsGrids.jl",
    format = Documenter.HTML(;
        canonical = "https://github.com/tensor4all/QuanticsGrids.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md", "API Reference" => "apireference.md"],
)

deploydocs(; repo = "github.com/tensor4all/QuanticsGrids.jl.git", devbranch = "main")
