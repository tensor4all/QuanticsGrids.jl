#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

using QuanticsGrids
using Chairmarks
using Random
using Printf

Random.seed!(42)

function create_equivalent_grids(::Val{D}, R, lower_bound, upper_bound; base=2, unfoldingscheme=:fused) where D
    old_grid = DiscretizedGrid{D}(R, lower_bound, upper_bound; base, unfoldingscheme)
    new_grid = NewDiscretizedGrid{D}(R; lower_bound, upper_bound, base, unfoldingscheme)
    return old_grid, new_grid
end

function generate_test_coordinates(grid, n_samples=1000)
    coord_samples = []
    D = typeof(grid) <: DiscretizedGrid ? length(grid.lower_bound) : ndims(grid)

    for _ in 1:n_samples
        coords = ntuple(D) do d
            lower = grid.lower_bound[d]
            upper = grid.upper_bound[d]
            lower + rand() * (upper - lower)
        end
        coords = D == 1 ? coords[1] : coords
        push!(coord_samples, coords)
    end
    return coord_samples
end

function generate_test_quantics(grid, n_samples=1000)
    quantics_samples = []

    for _ in 1:n_samples
        if typeof(grid) <: DiscretizedGrid
            D = length(grid.lower_bound)
            quantics = if grid.unfoldingscheme === :fused
                rand(1:grid.base^D, grid.R)
            else
                rand(1:grid.base, grid.R * D)
            end
        else
            local_dims = QuanticsGrids.localdimensions(grid)
            quantics = [rand(1:local_dim) for local_dim in local_dims]
        end
        push!(quantics_samples, quantics)
    end
    return quantics_samples
end

function benchmark_grid_config(config_name, D, R, lower_bound, upper_bound; base=2, unfoldingscheme=:fused, n_samples=1000)
    println("\n" * "="^60)
    println("$config_name")
    println("Dimensions: $D, R: $R, Base: $base, Unfolding: $unfoldingscheme")
    println("="^60)

    old_grid, new_grid = create_equivalent_grids(Val(D), R, lower_bound, upper_bound; base, unfoldingscheme)
    quantics_samples = generate_test_quantics(old_grid, n_samples)
    coord_samples = generate_test_coordinates(old_grid, n_samples)

    for i in 1:min(5, length(quantics_samples))
        quantics, coords = quantics_samples[i], coord_samples[i]

        old_result = quantics_to_origcoord(old_grid, quantics)
        new_result = quantics_to_origcoord(new_grid, quantics)

        if old_result isa Tuple && new_result isa Tuple
            all(isapprox.(old_result, new_result, rtol=1e-12)) || @warn "Results differ (quantics->origcoord)"
        elseif old_result isa Number && new_result isa Number
            isapprox(old_result, new_result, rtol=1e-12) || @warn "Results differ (quantics->origcoord)"
        end

        old_result_rev = origcoord_to_quantics(old_grid, coords)
        new_result_rev = origcoord_to_quantics(new_grid, coords)
        old_result_rev == new_result_rev || @warn "Results differ (origcoord->quantics)"
    end

    println("QUANTICS → ORIGCOORD")

    old_benchmark = @be begin
        for quantics in quantics_samples
            quantics_to_origcoord(old_grid, quantics)
        end
    end

    new_benchmark = @be begin
        for quantics in quantics_samples
            quantics_to_origcoord(new_grid, quantics)
        end
    end

    old_time = minimum([s.time for s in old_benchmark.samples]) / n_samples
    new_time = minimum([s.time for s in new_benchmark.samples]) / n_samples
    ratio = old_time / new_time

    @printf("  DiscretizedGrid:    %.1f ns/op\n", old_time * 1e9)
    @printf("  NewDiscretizedGrid: %.1f ns/op\n", new_time * 1e9)
    @printf("  Speedup: %.2fx %s\n", ratio, ratio > 1.1 ? "✅" : ratio > 0.9 ? "⚪" : "⚠️")

    println("\nORIGCOORD → QUANTICS")

    old_benchmark_rev = @be begin
        for coords in coord_samples
            origcoord_to_quantics(old_grid, coords)
        end
    end

    new_benchmark_rev = @be begin
        for coords in coord_samples
            origcoord_to_quantics(new_grid, coords)
        end
    end

    old_time_rev = minimum([s.time for s in old_benchmark_rev.samples]) / n_samples
    new_time_rev = minimum([s.time for s in new_benchmark_rev.samples]) / n_samples
    ratio_rev = old_time_rev / new_time_rev

    @printf("  DiscretizedGrid:    %.1f ns/op\n", old_time_rev * 1e9)
    @printf("  NewDiscretizedGrid: %.1f ns/op\n", new_time_rev * 1e9)
    @printf("  Speedup: %.2fx %s\n", ratio_rev, ratio_rev > 1.1 ? "✅" : ratio_rev > 0.9 ? "⚪" : "⚠️")

    return (
        old_time=old_time, new_time=new_time, ratio=ratio,
        old_time_rev=old_time_rev, new_time_rev=new_time_rev, ratio_rev=ratio_rev
    )
end
function main()
    println("QuanticsGrids Performance Benchmark")
    println("Julia $(VERSION), $(Threads.nthreads()) threads")

    configs = [
        ("1D Small", 1, 8, 0.0, 1.0, 2, :fused, 10000),
        ("1D Medium", 1, 16, -5.0, 5.0, 2, :fused, 5000),
        ("2D Fused", 2, 6, (0.0, 0.0), (1.0, 1.0), 2, :fused, 5000),
        ("2D Interleaved", 2, 6, (0.0, 0.0), (1.0, 1.0), 2, :interleaved, 5000),
        ("2D Medium", 2, 10, (-1.0, -1.0), (1.0, 1.0), 2, :fused, 2000),
        ("3D", 3, 5, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0), 2, :fused, 1000),
        ("Base-3 2D", 2, 6, (0.0, 0.0), (1.0, 1.0), 3, :fused, 2000),
        ("1D Large", 1, 20, 0.0, 1.0, 2, :fused, 1000)
    ]

    results = []
    for (name, D, R, lower, upper, base, scheme, n_samples) in configs
        result = benchmark_grid_config(name, D, R, lower, upper;
            base=base, unfoldingscheme=scheme, n_samples=n_samples)
        push!(results, (name, result))
    end

    print_summary(results)
    return results
end

function print_summary(results)
    println("\n" * "="^80)
    println("BENCHMARK SUMMARY")
    println("="^80)

    println("QUANTICS → ORIGCOORD")
    println(@sprintf("%-15s %12s %12s %12s", "Config", "Old (ns)", "New (ns)", "Speedup"))
    println("-"^80)

    for (name, result) in results
        @printf("%-15s %12.1f %12.1f %12.2fx\n",
            name, result.old_time * 1e9, result.new_time * 1e9, result.ratio)
    end

    println("\nORIGCOORD → QUANTICS")
    println(@sprintf("%-15s %12s %12s %12s", "Config", "Old (ns)", "New (ns)", "Speedup"))
    println("-"^80)

    for (name, result) in results
        @printf("%-15s %12.1f %12.1f %12.2fx\n",
            name, result.old_time_rev * 1e9, result.new_time_rev * 1e9, result.ratio_rev)
    end

    ratios = [r[2].ratio for r in results]
    ratios_rev = [r[2].ratio_rev for r in results]

    println("\n" * "="^80)
    println("PERFORMANCE ANALYSIS")
    println("="^80)

    @printf("Average speedup (quantics→origcoord): %.2fx\n", sum(ratios) / length(ratios))
    @printf("Average speedup (origcoord→quantics): %.2fx\n", sum(ratios_rev) / length(ratios_rev))

    faster_count = count(r -> r[2].ratio > 1.1, results)
    faster_count_rev = count(r -> r[2].ratio_rev > 1.1, results)

    println("\nNewDiscretizedGrid faster in $faster_count/$(length(results)) forward tests")
    println("NewDiscretizedGrid faster in $faster_count_rev/$(length(results)) reverse tests")

    overall_improvement = (sum(ratios) + sum(ratios_rev)) / (2 * length(results))
    if overall_improvement > 1.1
        println("\n🎉 NewDiscretizedGrid shows significant performance improvements!")
    elseif overall_improvement > 0.9
        println("\n✅ NewDiscretizedGrid maintains equivalent performance")
    else
        println("\n⚠️  NewDiscretizedGrid shows performance regression")
    end
end

main()
