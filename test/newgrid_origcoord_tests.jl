@testitem "origcoord functions - basic 1D grid" begin
    grid = QuanticsGrids.NewDiscretizedGrid((4,); lower_bound=(0.0,), upper_bound=(1.0,))

    # Test grididx <-> origcoord conversion
    @test QuanticsGrids.grididx_to_origcoord(grid, (1,)) == 0.0
    @test QuanticsGrids.grididx_to_origcoord(grid, (2^4,)) ≈ 1.0 - 1.0 / 2^4

    # Test origcoord <-> grididx conversion  
    @test QuanticsGrids.origcoord_to_grididx(grid, 0.0) == 1
    @test QuanticsGrids.origcoord_to_grididx(grid, 0.5) == 2^3 + 1 # 0b1000 + 1

    # Test round-trip conversions
    for grididx in [(1,), (5,), (16,)]
        origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
        @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == only(grididx)
    end

    # Test quantics <-> origcoord conversions (should work through existing functions)
    quantics = [1, 1, 1, 2]  # should correspond to grididx = 2
    expected_origcoord = 1.0 / 2^4
    @test QuanticsGrids.quantics_to_origcoord(grid, quantics) ≈ expected_origcoord
    @test QuanticsGrids.origcoord_to_quantics(grid, expected_origcoord) == quantics
end

@testitem "origcoord functions - 2D grid" begin
    grid = QuanticsGrids.NewDiscretizedGrid((3, 4); lower_bound=(0.0, 1.0), upper_bound=(2.0, 3.0))

    # Test boundary values
    @test QuanticsGrids.grididx_to_origcoord(grid, (1, 1)) == (0.0, 1.0)
    @test all(isapprox.(QuanticsGrids.grididx_to_origcoord(grid, (2^3, 2^4)), (2.0 - 2.0 / 2^3, 3.0 - 2.0 / 2^4)))

    # Test middle values
    mid_grididx = (2^2 + 1, 2^3 + 1)  # Middle of each dimension
    expected_origcoord = (0.0 + 2.0 * 2^2 / 2^3, 1.0 + 2.0 * 2^3 / 2^4)
    @test all(isapprox.(QuanticsGrids.grididx_to_origcoord(grid, mid_grididx), expected_origcoord))

    # Test round-trip conversions
    for _ in 1:20
        grididx = (rand(1:2^3), rand(1:2^4))
        origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
        @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx
    end
end

@testitem "origcoord functions - different base" begin
    base = 3
    grid = QuanticsGrids.NewDiscretizedGrid((2, 3); base=base, lower_bound=(-1.0, 0.0), upper_bound=(1.0, 6.0))

    # Test boundary values  
    @test QuanticsGrids.grididx_to_origcoord(grid, (1, 1)) == (-1.0, 0.0)
    @test all(isapprox.(QuanticsGrids.grididx_to_origcoord(grid, (base^2, base^3)), (1.0 - 2.0 / base^2, 6.0 - 6.0 / base^3)))

    # Test round-trip conversions
    for _ in 1:20
        grididx = (rand(1:base^2), rand(1:base^3))
        origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
        @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx
    end

    # Test quantics <-> origcoord conversions
    for _ in 1:20
        quantics = rand(1:base, length(grid))
        grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
        origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)

        @test all(isapprox.(QuanticsGrids.quantics_to_origcoord(grid, quantics), origcoord))
        @test QuanticsGrids.origcoord_to_quantics(grid, origcoord) == quantics
    end
end

@testitem "origcoord functions - boundary checking" begin
    grid = QuanticsGrids.NewDiscretizedGrid((2, 2); lower_bound=(0.0, 0.0), upper_bound=(1.0, 1.0))

    # Test coordinates within bounds
    @test QuanticsGrids.origcoord_to_grididx(grid, (0.0, 0.0)) == (1, 1)
    @test QuanticsGrids.origcoord_to_grididx(grid, (0.5, 0.5)) == (3, 3)
    @test QuanticsGrids.origcoord_to_grididx(grid, (0.999, 0.999)) == (4, 4)

    # Test coordinates outside bounds should throw errors
    @test_throws BoundsError QuanticsGrids.origcoord_to_grididx(grid, (-0.1, 0.5))
    @test_throws BoundsError QuanticsGrids.origcoord_to_grididx(grid, (0.5, -0.1))
    @test_throws BoundsError QuanticsGrids.origcoord_to_grididx(grid, (1.1, 0.5))
    @test_throws BoundsError QuanticsGrids.origcoord_to_grididx(grid, (0.5, 1.1))
    @test_throws BoundsError QuanticsGrids.origcoord_to_grididx(grid, (-0.1, -0.1))
    @test_throws BoundsError QuanticsGrids.origcoord_to_grididx(grid, (1.1, 1.1))
end

@testitem "origcoord functions - fused indices grid" begin
    # Test with a more complex grid that has fused indices
    grid = QuanticsGrids.NewDiscretizedGrid(
        (:x, :y),
        [[(:x, 2), (:y, 1)], [(:x, 1)], [(:y, 2)]];
        lower_bound=(0.0, -1.0),
        upper_bound=(4.0, 1.0)
    )

    # Test some conversions
    for _ in 1:50
        quantics = [rand(1:4), rand(1:2), rand(1:2)]  # Based on sitedim values
        grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
        origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)

        # Test round-trip conversion
        @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx
        @test all(isapprox.(QuanticsGrids.quantics_to_origcoord(grid, quantics), origcoord))
        @test QuanticsGrids.origcoord_to_quantics(grid, origcoord) == quantics
    end
end

@testitem "origcoord functions - stress test" begin
    # Test high-dimensional grids
    grid = QuanticsGrids.NewDiscretizedGrid((2, 3, 2, 4); base=2,
        lower_bound=(0.0, -2.0, 1.0, -5.0),
        upper_bound=(1.0, 2.0, 3.0, 5.0))

    for _ in 1:30
        # Test random grid indices
        grididx = (rand(1:2^2), rand(1:2^3), rand(1:2^2), rand(1:2^4))
        origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
        @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx

        # Test random quantics
        quantics = rand(1:2, length(grid))
        grididx2 = QuanticsGrids.quantics_to_grididx(grid, quantics)
        origcoord2 = QuanticsGrids.grididx_to_origcoord(grid, grididx2)

        @test all(isapprox.(QuanticsGrids.quantics_to_origcoord(grid, quantics), origcoord2))
        @test QuanticsGrids.origcoord_to_quantics(grid, origcoord2) == quantics
    end
end

@testitem "edge cases - numerical precision and floating point" begin
    # Test with very small grid steps that might cause floating point issues
    grid = QuanticsGrids.NewDiscretizedGrid((30, 25);
        lower_bound=(1e-15, -1e-15),
        upper_bound=(1e-14, 1e-14))

    # Test conversions near boundaries where floating point precision matters
    for _ in 1:20
        grididx = (rand(1:2^30), rand(1:2^25))
        origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
        recovered_grididx = QuanticsGrids.origcoord_to_grididx(grid, origcoord)
        @test recovered_grididx == grididx
    end

    # Test coordinates very close to boundaries
    eps_coord_x = 1e-15 + 1e-20  # Just above lower bound
    eps_coord_y = 1e-14 - 1e-20  # Just below upper bound
    @test QuanticsGrids.origcoord_to_grididx(grid, (eps_coord_x, eps_coord_y)) isa Tuple
end

@testitem "edge cases - maximum R values and large numbers" begin
    # Test with largest R values that still fit in computation
    # R=60 gives 2^60 ≈ 1.15e18, near but below Int64 max (9.22e18)
    large_R = 58  # Be conservative to avoid overflow
    grid = QuanticsGrids.NewDiscretizedGrid((large_R,))

    # Test extreme indices
    min_quantics = ones(Int, large_R)
    max_quantics = fill(2, large_R)

    @test QuanticsGrids.quantics_to_grididx(grid, min_quantics) == 1
    @test QuanticsGrids.quantics_to_grididx(grid, max_quantics) == 2^large_R
    @test QuanticsGrids.grididx_to_quantics(grid, (1,)) == min_quantics
    @test QuanticsGrids.grididx_to_quantics(grid, (2^large_R,)) == max_quantics

    # Test some random middle values
    for _ in 1:10
        quantics = rand(1:2, large_R)
        grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
        @test QuanticsGrids.grididx_to_quantics(grid, grididx) == quantics
    end
end

@testitem "edge cases - degenerate grids with R=0" begin
    # Grid with all R=0 (should only have one grid point each)
    grid = QuanticsGrids.NewDiscretizedGrid((0, 0, 0))

    # Only valid grid index should be (1, 1, 1)
    @test QuanticsGrids.grididx_to_quantics(grid, (1, 1, 1)) == Int[]
    @test QuanticsGrids.quantics_to_grididx(grid, Int[]) == (1, 1, 1)

    # Test origcoord functions
    @test QuanticsGrids.grididx_to_origcoord(grid, (1, 1, 1)) == (0.0, 0.0, 0.0)
    @test QuanticsGrids.origcoord_to_grididx(grid, (0.0, 0.0, 0.0)) == (1, 1, 1)
end

@testitem "edge cases - mixed R values including zeros" begin
    # Grid with mix of zero and non-zero R values
    grid = QuanticsGrids.NewDiscretizedGrid((0, 5, 0, 3, 0))

    # Test conversion with valid quantics length
    quantics = [1, 2, 1, 2, 1, 2, 1, 2]  # Only bits for dimensions with R>0
    grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
    @test grididx[1] == 1  # R=0 dimensions should always be 1
    @test grididx[3] == 1  # R=0 dimensions should always be 1
    @test grididx[5] == 1  # R=0 dimensions should always be 1
    @test QuanticsGrids.grididx_to_quantics(grid, grididx) == quantics

    # Test boundary conditions
    for _ in 1:20
        valid_quantics = rand(1:2, length(grid))
        grididx_local = QuanticsGrids.quantics_to_grididx(grid, valid_quantics)
        @test grididx_local[1] == 1
        @test grididx_local[3] == 1
        @test grididx_local[5] == 1
        @test 1 <= grididx_local[2] <= 2^5
        @test 1 <= grididx_local[4] <= 2^3
    end
end

@testitem "edge cases - extreme base values" begin
    # Test with maximum reasonable base
    max_base = 16
    grid = QuanticsGrids.NewDiscretizedGrid((3, 4); base=max_base)

    for _ in 1:20
        quantics = rand(1:max_base, length(grid))
        grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
        @test QuanticsGrids.grididx_to_quantics(grid, grididx) == quantics
        @test all(1 .<= grididx .<= (max_base .^ grid.Rs))
    end

    # Test with base=2 (edge case, minimum valid base)
    grid2 = QuanticsGrids.NewDiscretizedGrid((10, 8); base=2)
    for _ in 1:20
        quantics = rand(1:2, length(grid2))
        grididx = QuanticsGrids.quantics_to_grididx(grid2, quantics)
        @test QuanticsGrids.grididx_to_quantics(grid2, grididx) == quantics
    end
end

@testitem "edge cases - extreme coordinate ranges" begin
    # Test with very large coordinate ranges
    grid = QuanticsGrids.NewDiscretizedGrid((10, 8);
        lower_bound=(-1e10, -1e15),
        upper_bound=(1e10, 1e15))

    for _ in 1:20
        grididx = (rand(1:2^10), rand(1:2^8))
        origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
        @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx
    end

    # Test with very small coordinate ranges
    grid2 = QuanticsGrids.NewDiscretizedGrid((5, 6);
        lower_bound=(1e-12, -1e-12),
        upper_bound=(1e-12 + 1e-15, -1e-12 + 1e-15))

    for _ in 1:10
        grididx = (rand(1:2^5), rand(1:2^6))
        origcoord = QuanticsGrids.grididx_to_origcoord(grid2, grididx)
        recovered = QuanticsGrids.origcoord_to_grididx(grid2, origcoord)
        @test recovered == grididx
    end
end

@testitem "edge cases - asymmetric ranges with negative bounds" begin
    # Test grids with very asymmetric and negative bounds
    grid = QuanticsGrids.NewDiscretizedGrid((4, 6, 3);
        lower_bound=(-1000.0, -0.001, -1e6),
        upper_bound=(-999.0, 0.001, 1e6))

    for _ in 1:30
        quantics = rand(1:2, length(grid))
        grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
        origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)

        # Verify coordinates are within bounds
        @test all(QuanticsGrids.grid_min(grid) .<= origcoord .<= QuanticsGrids.grid_max(grid))

        # Test round-trip
        @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx
        @test all(isapprox.(QuanticsGrids.quantics_to_origcoord(grid, quantics), origcoord))
    end
end

@testitem "stress test - maximum complexity fused indices" begin
    # Create maximally complex fused index structure
    # All variables distributed unevenly across sites
    grid = QuanticsGrids.NewDiscretizedGrid(
        (:a, :b, :c, :d, :e, :f, :g, :h),
        [
            [(:a, 8), (:b, 7), (:c, 6), (:d, 5)],     # 4 vars, site dim = base^4
            [(:e, 4), (:f, 3)],                        # 2 vars, site dim = base^2  
            [(:g, 2)],                                  # 1 var
            [(:h, 1)],                                  # 1 var
            [(:a, 7), (:c, 5), (:e, 3), (:g, 1)],     # 4 vars, complex pattern
            [(:b, 6), (:d, 4), (:f, 2)],              # 3 vars
            [(:a, 6), (:h, 2)],                        # 2 vars from distant indices
            [(:b, 5), (:c, 4), (:d, 3)],              # 3 vars
            [(:e, 2), (:f, 1)],                        # 2 vars
            [(:a, 5), (:g, 3)],                        # 2 vars, distant indices
            [(:b, 4), (:c, 3), (:d, 2), (:h, 3)],     # 4 vars
            [(:a, 4), (:e, 1)],                        # 2 vars
            [(:b, 3), (:f, 4)],                        # 2 vars
            [(:a, 3), (:c, 2)],                        # 2 vars
            [(:b, 2), (:d, 1)],                        # 2 vars
            [(:a, 2)],                                  # 1 var
            [(:b, 1)],                                  # 1 var
            [(:a, 1)],                                  # 1 var
            [(:c, 1)]                                   # 1 var
        ];
        base=2
    )

    # Stress test with many conversions
    for _ in 1:100
        # Generate valid quantics based on site dimensions
        quantics = [
            rand(1:16),   # site 1: 2^4
            rand(1:4),    # site 2: 2^2
            rand(1:2),    # site 3: 2^1
            rand(1:2),    # site 4: 2^1
            rand(1:16),   # site 5: 2^4
            rand(1:8),    # site 6: 2^3
            rand(1:4),    # site 7: 2^2
            rand(1:8),    # site 8: 2^3
            rand(1:4),    # site 9: 2^2
            rand(1:4),    # site 10: 2^2
            rand(1:16),   # site 11: 2^4
            rand(1:4),    # site 12: 2^2
            rand(1:4),    # site 13: 2^2
            rand(1:4),    # site 14: 2^2
            rand(1:4),    # site 15: 2^2
            rand(1:2),    # site 16: 2^1
            rand(1:2),    # site 17: 2^1
            rand(1:2),    # site 18: 2^1
            rand(1:2)     # site 19: 2^1
        ]

        grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
        recovered = QuanticsGrids.grididx_to_quantics(grid, grididx)
        @test recovered == quantics

        # Test origcoord functions too
        origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
        @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx
    end
end

@testitem "boundary stress test - coordinates exactly on boundaries" begin
    grid = QuanticsGrids.NewDiscretizedGrid((8, 6);
        lower_bound=(-2.0, 3.0),
        upper_bound=(5.0, 8.0))

    # Test coordinates exactly at lower bounds
    @test QuanticsGrids.origcoord_to_grididx(grid, (-2.0, 3.0)) == (1, 1)

    # Test coordinates exactly at effective upper bounds (grid_max)
    grid_max_coord = (-2.0 + 7.0 * (2^8 - 1) / 2^8, 3.0 + 5.0 * (2^6 - 1) / 2^6)
    max_grididx = (2^8, 2^6)
    @test QuanticsGrids.origcoord_to_grididx(grid, grid_max_coord) == max_grididx

    # Test coordinates just inside boundaries
    eps = 1e-12
    inside_lower = (-2.0 + eps, 3.0 + eps)
    inside_upper = (grid_max_coord[1] - eps, grid_max_coord[2] - eps)

    @test QuanticsGrids.origcoord_to_grididx(grid, inside_lower) isa Tuple
    @test QuanticsGrids.origcoord_to_grididx(grid, inside_upper) isa Tuple

    # Test coordinates just outside boundaries should error
    outside_lower = (-2.0 - eps, 3.0 - eps)
    outside_upper = (5.0 + eps, 8.0 + eps)

    @test_throws BoundsError QuanticsGrids.origcoord_to_grididx(grid, outside_lower)
    @test_throws BoundsError QuanticsGrids.origcoord_to_grididx(grid, outside_upper)
end

@testitem "performance stress test - rapid conversions" begin
    # Test performance with moderately sized grids and many conversions
    grid = QuanticsGrids.NewDiscretizedGrid((12, 10, 8, 6); base=3)

    # Pre-generate test data to avoid timing allocation
    n_tests = 1000
    test_quantics = [rand(1:3, length(grid)) for _ in 1:n_tests]
    test_grididx = [(rand(1:3^12), rand(1:3^10), rand(1:3^8), rand(1:3^6)) for _ in 1:n_tests]

    # Test quantics -> grididx -> quantics round trips
    for i in 1:n_tests
        quantics = test_quantics[i]
        grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
        @test QuanticsGrids.grididx_to_quantics(grid, grididx) == quantics
    end

    # Test grididx -> quantics -> grididx round trips  
    for i in 1:n_tests
        grididx = test_grididx[i]
        quantics = QuanticsGrids.grididx_to_quantics(grid, grididx)
        @test QuanticsGrids.quantics_to_grididx(grid, quantics) == grididx
    end

    # Test origcoord conversions
    for i in 1:min(n_tests, 200)  # Fewer origcoord tests since they're slower
        grididx = test_grididx[i]
        origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
        @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx
    end
end
