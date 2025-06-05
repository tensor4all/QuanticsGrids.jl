# Benchmarking quantics_to_grididx and grididx_to_quantics for NewDiscretizedGrid
# Usage: julia --project -e 'using Pkg; Pkg.add("Chairmarks"); include("bench_newgrid_quantics.jl")'

using QuanticsGrids
using Random
using Chairmarks

function random_quantics(grid)
    [rand(1:QuanticsGrids.sitedim(grid, site)) for site in eachindex(grid.indextable)]
end

function random_grididx(grid)
    ntuple(d -> rand(1:grid.base^grid.Rs[d]), ndims(grid))
end

function make_grids()
    [
        QuanticsGrids.NewDiscretizedGrid((20,)),
        QuanticsGrids.NewDiscretizedGrid((10, 12)),
        QuanticsGrids.NewDiscretizedGrid((6, 7, 8)),
        QuanticsGrids.NewDiscretizedGrid((:a, :b, :c), [[(:a, 4)], [(:a, 3)], [(:a, 2)], [(:a, 1)], [(:b, 1)], [(:b, 2)], [(:b, 3)], [(:c, 1)], [(:c, 2)], [(:c, 3)]]),
        QuanticsGrids.NewDiscretizedGrid((5, 3, 17); base=13),
        QuanticsGrids.NewDiscretizedGrid(ntuple(i -> 4 + (i % 3), 8); base=3)
    ]
end

grids = make_grids()

println("Running Chairmarks benchmarks for NewDiscretizedGrid quantics/grididx conversions...")

for (i, grid) in enumerate(grids)
    quantics = random_quantics(grid)
    grididx = random_grididx(grid)

    println("Grid $(i): $(typeof(grid).name.name)")
    println("  quantics_to_grididx: ", @b QuanticsGrids.quantics_to_grididx($grid, $quantics))
    println("  grididx_to_quantics: ", @b QuanticsGrids.grididx_to_quantics($grid, $grididx))
    println()
end
