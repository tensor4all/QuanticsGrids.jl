#!/usr/bin/env julia
"""
Benchmark script comparing quantics_to_origcoord performance between 
NewDiscretizedGrid and DiscretizedGrid implementations.

This script tests various grid configurations to ensure the new implementation
doesn't compromise on performance.
"""

using Pkg
Pkg.activate(@__DIR__)

using QuanticsGrids
using Chairmarks
using Random
using Printf

# Set random seed for reproducibility
Random.seed!(42)

"""
Create equivalent DiscretizedGrid and NewDiscretizedGrid for benchmarking
"""
function create_equivalent_grids(::Val{D}, R, lower_bound, upper_bound; base=2, unfoldingscheme=:fused) where D
    # Create DiscretizedGrid
    old_grid = DiscretizedGrid{D}(R, lower_bound, upper_bound; base, unfoldingscheme)

    # Create NewDiscretizedGrid with equivalent parameters
    new_grid = NewDiscretizedGrid{D}(R; lower_bound, upper_bound, base, unfoldingscheme)

    return old_grid, new_grid
end

"""
Generate random coordinate vectors for reverse direction testing
"""
function generate_test_coordinates(grid, n_samples=1000)
    coord_samples = []

    for _ in 1:n_samples
        if typeof(grid) <: DiscretizedGrid
            D = length(grid.lower_bound)
            # Generate random coordinates within bounds
            coords = ntuple(D) do d
                lower = grid.lower_bound[d]
                upper = grid.upper_bound[d]
                lower + rand() * (upper - lower)
            end
            if D == 1
                coords = coords[1]  # Convert to scalar for 1D
            end
        else  # NewDiscretizedGrid
            D = ndims(grid)
            coords = ntuple(D) do d
                lower = grid.lower_bound[d]
                upper = grid.upper_bound[d]
                lower + rand() * (upper - lower)
            end
            if D == 1
                coords = coords[1]  # Convert to scalar for 1D
            end
        end
        push!(coord_samples, coords)
    end

    return coord_samples
end

"""
Generate random quantics vectors for testing
"""
function generate_test_quantics(grid, n_samples=1000)
    quantics_samples = []

    for _ in 1:n_samples
        if typeof(grid) <: DiscretizedGrid
            # For old grid, use the global dimension structure
            D = length(grid.lower_bound)  # Get dimension from the lower_bound field
            if grid.unfoldingscheme === :fused
                quantics = rand(1:grid.base^D, grid.R)
            else  # interleaved
                quantics = rand(1:grid.base, grid.R * D)
            end
        else  # NewDiscretizedGrid
            # For new grid, use the localdimensions
            local_dims = QuanticsGrids.localdimensions(grid)
            quantics = [rand(1:local_dim) for local_dim in local_dims]
        end
        push!(quantics_samples, quantics)
    end

    return quantics_samples
end

"""
Run benchmark for a specific grid configuration
"""
function benchmark_grid_config(config_name, D, R, lower_bound, upper_bound; base=2, unfoldingscheme=:fused, n_samples=1000)
    println("\n" * "="^60)
    println("Benchmarking: $config_name")
    println("Dimensions: $D, R: $R, Base: $base, Unfolding: $unfoldingscheme")
    println("Domain: $lower_bound → $upper_bound")
    println("="^60)

    # Create grids
    old_grid, new_grid = create_equivalent_grids(Val(D), R, lower_bound, upper_bound; base, unfoldingscheme)

    # Generate test data - use old grid structure for compatibility
    quantics_samples = generate_test_quantics(old_grid, n_samples)
    coord_samples = generate_test_coordinates(old_grid, n_samples)

    println("Grid Info:")
    println("  Old grid total points: $(base^R)^$D = $(base^(R*D))")
    if typeof(new_grid) <: NewDiscretizedGrid
        new_total = prod(base .^ new_grid.Rs)
        println("  New grid total points: $new_total")
    end
    println("  Number of test samples: $n_samples")

    # Verify equivalence for a few samples
    println("\nVerifying equivalence...")
    n_verify = min(10, length(quantics_samples))
    for i in 1:n_verify
        quantics = quantics_samples[i]
        coords = coord_samples[i]
        try
            # Test quantics -> origcoord
            old_result = quantics_to_origcoord(old_grid, quantics)
            new_result = quantics_to_origcoord(new_grid, quantics)

            if old_result isa Tuple && new_result isa Tuple
                if !all(isapprox.(old_result, new_result, rtol=1e-12))
                    @warn "quantics->origcoord results differ for sample $i: old=$old_result, new=$new_result"
                end
            elseif old_result isa Number && new_result isa Number
                if !isapprox(old_result, new_result, rtol=1e-12)
                    @warn "quantics->origcoord results differ for sample $i: old=$old_result, new=$new_result"
                end
            else
                @warn "quantics->origcoord result types differ for sample $i: old=$(typeof(old_result)), new=$(typeof(new_result))"
            end

            # Test origcoord -> quantics
            old_result_rev = origcoord_to_quantics(old_grid, coords)
            new_result_rev = origcoord_to_quantics(new_grid, coords)

            if old_result_rev != new_result_rev
                @show coords
                @warn "origcoord->quantics results differ for sample $i: old=$old_result_rev, new=$new_result_rev"
            end
        catch e
            @warn "Error verifying sample $i: $e"
        end
    end

    # Benchmark quantics -> origcoord direction
    println("\n" * "="^40)
    println("QUANTICS → ORIGCOORD")
    println("="^40)

    # Benchmark old implementation
    println("Benchmarking DiscretizedGrid...")
    old_benchmark = @be begin
        for quantics in $quantics_samples
            quantics_to_origcoord($old_grid, quantics)
        end
    end

    # Benchmark new implementation  
    println("Benchmarking NewDiscretizedGrid...")
    new_benchmark = @be begin
        for quantics in $quantics_samples
            quantics_to_origcoord($new_grid, quantics)
        end
    end

    # Calculate performance ratio (time per conversion)
    old_time_total = minimum([s.time for s in old_benchmark.samples])
    new_time_total = minimum([s.time for s in new_benchmark.samples])
    old_time = old_time_total / n_samples  # Time per individual conversion
    new_time = new_time_total / n_samples  # Time per individual conversion
    ratio = old_time / new_time

    # Print results (summary only, without detailed samples)
    println("\nResults:")
    @printf("  DiscretizedGrid     (old): %.3f ns per conversion\n", old_time * 1e9)
    @printf("  NewDiscretizedGrid  (new): %.3f ns per conversion\n", new_time * 1e9)
    @printf("  Performance ratio (old/new): %.3fx\n", ratio)

    if ratio > 1.1
        println("  ✅ New implementation is faster")
    elseif ratio > 0.9
        println("  ⚪ Performance is similar")
    else
        println("  ⚠️  New implementation is slower")
    end

    # Benchmark origcoord -> quantics direction
    println("\n" * "="^40)
    println("ORIGCOORD → QUANTICS")
    println("="^40)

    # Benchmark old implementation (reverse)
    println("Benchmarking DiscretizedGrid...")
    old_benchmark_rev = @be begin
        for coords in $coord_samples
            origcoord_to_quantics($old_grid, coords)
        end
    end

    # Benchmark new implementation (reverse)
    println("Benchmarking NewDiscretizedGrid...")
    new_benchmark_rev = @be begin
        for coords in $coord_samples
            origcoord_to_quantics($new_grid, coords)
        end
    end

    # Calculate performance ratio for reverse direction
    old_time_total_rev = minimum([s.time for s in old_benchmark_rev.samples])
    new_time_total_rev = minimum([s.time for s in new_benchmark_rev.samples])
    old_time_rev = old_time_total_rev / n_samples
    new_time_rev = new_time_total_rev / n_samples
    ratio_rev = old_time_rev / new_time_rev

    # Print reverse results
    println("\nResults:")
    @printf("  DiscretizedGrid     (old): %.3f ns per conversion\n", old_time_rev * 1e9)
    @printf("  NewDiscretizedGrid  (new): %.3f ns per conversion\n", new_time_rev * 1e9)
    @printf("  Performance ratio (old/new): %.3fx\n", ratio_rev)

    if ratio_rev > 1.1
        println("  ✅ New implementation is faster")
    elseif ratio_rev > 0.9
        println("  ⚪ Performance is similar")
    else
        println("  ⚠️  New implementation is slower")
    end

    return (
        old_time=old_time, new_time=new_time, ratio=ratio,
        old_time_rev=old_time_rev, new_time_rev=new_time_rev, ratio_rev=ratio_rev
    )
end

"""
Main benchmarking suite
"""
function main()
    println("QuanticsGrids Benchmark: quantics_to_origcoord performance comparison")
    println("Julia Version: $(VERSION)")
    println("Threads: $(Threads.nthreads())")

    results = []

    # Test configuration 1: Small 1D grid
    result1 = benchmark_grid_config(
        "Small 1D Grid", 1, 8, 0.0, 1.0;
        base=2, unfoldingscheme=:fused, n_samples=10000
    )
    push!(results, ("Small 1D", result1))

    # Test configuration 2: Medium 1D grid
    result2 = benchmark_grid_config(
        "Medium 1D Grid", 1, 16, -5.0, 5.0;
        base=2, unfoldingscheme=:fused, n_samples=5000
    )
    push!(results, ("Medium 1D", result2))

    # Test configuration 3: Small 2D grid (fused)
    result3 = benchmark_grid_config(
        "Small 2D Grid (Fused)", 2, 6, (0.0, 0.0), (1.0, 1.0);
        base=2, unfoldingscheme=:fused, n_samples=5000
    )
    push!(results, ("Small 2D Fused", result3))

    # Test configuration 4: Small 2D grid (interleaved)
    result4 = benchmark_grid_config(
        "Small 2D Grid (Interleaved)", 2, 6, (0.0, 0.0), (1.0, 1.0);
        base=2, unfoldingscheme=:interleaved, n_samples=5000
    )
    push!(results, ("Small 2D Interleaved", result4))

    # Test configuration 5: Medium 2D grid
    result5 = benchmark_grid_config(
        "Medium 2D Grid", 2, 10, (-1.0, -1.0), (1.0, 1.0);
        base=2, unfoldingscheme=:fused, n_samples=2000
    )
    push!(results, ("Medium 2D", result5))

    # Test configuration 6: 3D grid
    result6 = benchmark_grid_config(
        "3D Grid", 3, 5, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0);
        base=2, unfoldingscheme=:fused, n_samples=1000
    )
    push!(results, ("3D", result6))

    # Test configuration 7: Base-3 grid
    result7 = benchmark_grid_config(
        "Base-3 2D Grid", 2, 6, (0.0, 0.0), (1.0, 1.0);
        base=3, unfoldingscheme=:fused, n_samples=2000
    )
    push!(results, ("Base-3 2D", result7))

    # Test configuration 8: Large 1D grid (stress test)
    result8 = benchmark_grid_config(
        "Large 1D Grid", 1, 20, 0.0, 1.0;
        base=2, unfoldingscheme=:fused, n_samples=1000
    )
    push!(results, ("Large 1D", result8))

    # Summary
    println("\n" * "="^100)
    println("BENCHMARK SUMMARY")
    println("="^100)
    println("QUANTICS → ORIGCOORD")
    println(@sprintf("%-20s %15s %15s %15s", "Configuration", "Old Time (ns)", "New Time (ns)", "Ratio (old/new)"))
    println("-"^100)

    for (name, result) in results
        @printf("%-20s %15.3f %15.3f %15.3fx\n",
            name, result.old_time * 1e9, result.new_time * 1e9, result.ratio)
    end

    # Overall statistics for forward direction
    ratios = [r[2].ratio for r in results]
    avg_ratio = sum(ratios) / length(ratios)
    min_ratio = minimum(ratios)
    max_ratio = maximum(ratios)

    println("-"^100)
    @printf("%-20s %15s %15s %15.3fx\n", "Average", "", "", avg_ratio)
    @printf("%-20s %15s %15s %15.3fx\n", "Min", "", "", min_ratio)
    @printf("%-20s %15s %15s %15.3fx\n", "Max", "", "", max_ratio)

    # Summary for reverse direction
    println("\n" * "="^100)
    println("ORIGCOORD → QUANTICS")
    println(@sprintf("%-20s %15s %15s %15s", "Configuration", "Old Time (ns)", "New Time (ns)", "Ratio (old/new)"))
    println("-"^100)

    for (name, result) in results
        @printf("%-20s %15.3f %15.3f %15.3fx\n",
            name, result.old_time_rev * 1e9, result.new_time_rev * 1e9, result.ratio_rev)
    end

    # Overall statistics for reverse direction
    ratios_rev = [r[2].ratio_rev for r in results]
    avg_ratio_rev = sum(ratios_rev) / length(ratios_rev)
    min_ratio_rev = minimum(ratios_rev)
    max_ratio_rev = maximum(ratios_rev)

    println("-"^100)
    @printf("%-20s %15s %15s %15.3fx\n", "Average", "", "", avg_ratio_rev)
    @printf("%-20s %15s %15s %15.3fx\n", "Min", "", "", min_ratio_rev)
    @printf("%-20s %15s %15s %15.3fx\n", "Max", "", "", max_ratio_rev)

    println("\nPerformance Analysis:")

    # Forward direction analysis
    faster_count = count(r -> r[2].ratio > 1.1, results)
    similar_count = count(r -> 0.9 <= r[2].ratio <= 1.1, results)
    slower_count = count(r -> r[2].ratio < 0.9, results)

    println("QUANTICS → ORIGCOORD:")
    println("  ✅ New implementation faster: $faster_count/$(length(results))")
    println("  ⚪ Performance similar:       $similar_count/$(length(results))")
    println("  ⚠️  New implementation slower: $slower_count/$(length(results))")

    # Reverse direction analysis
    faster_count_rev = count(r -> r[2].ratio_rev > 1.1, results)
    similar_count_rev = count(r -> 0.9 <= r[2].ratio_rev <= 1.1, results)
    slower_count_rev = count(r -> r[2].ratio_rev < 0.9, results)

    println("\nORIGCOORD → QUANTICS:")
    println("  ✅ New implementation faster: $faster_count_rev/$(length(results))")
    println("  ⚪ Performance similar:       $similar_count_rev/$(length(results))")
    println("  ⚠️  New implementation slower: $slower_count_rev/$(length(results))")

    # Combined analysis
    avg_combined = (avg_ratio + avg_ratio_rev) / 2

    if avg_ratio >= 1.1
        println("\n🎉 Overall: NewDiscretizedGrid shows better performance!")
    elseif avg_ratio >= 0.9
        println("\n✅ Overall: NewDiscretizedGrid maintains equivalent performance!")
    else
        println("\n⚠️  Overall: NewDiscretizedGrid shows worse performance.")
    end

    return results
end

# Run the benchmarks if this script is executed directly
# if abspath(PROGRAM_FILE) == @__FILE__
main()
# end
