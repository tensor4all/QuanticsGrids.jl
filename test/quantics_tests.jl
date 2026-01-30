@testitem "constructor, square grid" begin
    grid = DiscretizedGrid{2}(10; unfoldingscheme=:interleaved)
    @test QuanticsGrids.grid_Rs(grid) == (10, 10)
    @test only(QuanticsGrids.grid_indextable(grid)[1])[1] == Symbol(1)
    @test only(QuanticsGrids.grid_indextable(grid)[1])[2] == 1
    @test only(QuanticsGrids.grid_indextable(grid)[2])[1] == Symbol(2)
    @test only(QuanticsGrids.grid_indextable(grid)[2])[2] == 1
end

@testitem "constructor, rectangular grid" begin
    grid = DiscretizedGrid((3, 5); unfoldingscheme=:interleaved)
    @test only(QuanticsGrids.grid_indextable(grid)[6])[1] == Symbol(2)
    @test only(QuanticsGrids.grid_indextable(grid)[6])[2] == 3
    @test only(QuanticsGrids.grid_indextable(grid)[7])[1] == Symbol(2)
    @test only(QuanticsGrids.grid_indextable(grid)[7])[2] == 4
    @test only(QuanticsGrids.grid_indextable(grid)[8])[1] == Symbol(2)
    @test only(QuanticsGrids.grid_indextable(grid)[8])[2] == 5
end

@testitem "quantics_to_grididx, rectangular grid" begin
    grid = DiscretizedGrid((3, 5); unfoldingscheme=:interleaved)
    @test quantics_to_grididx(grid, [1, 2, 1, 2, 1, 2, 1, 2]) == (1, 30)
end

@testitem "grididx_to_quantics, rectangular grid" begin
    grid = DiscretizedGrid((3, 5); unfoldingscheme=:interleaved)
    @test grididx_to_quantics(grid, (1, 30)) == [1, 2, 1, 2, 1, 2, 1, 2]
end


@testitem "grouped unfoldingscheme quantics conversions" begin
    Rs = (2, 1)
    grid = DiscretizedGrid(Rs; unfoldingscheme=:grouped)
    @test QuanticsGrids.localdimensions(grid) == fill(2, sum(Rs))

    grididx = (3, 2)
    expected_quantics = [2, 1, 2]
    @test grididx_to_quantics(grid, grididx) == expected_quantics
    @test quantics_to_grididx(grid, expected_quantics) == grididx

    grid_base3 = DiscretizedGrid(Rs; base=3, unfoldingscheme=:grouped)
    @test QuanticsGrids.localdimensions(grid_base3) == fill(3, sum(Rs))

    grididx_base3 = (5, 3)
    expected_quantics_base3 = [2, 2, 3]
    @test grididx_to_quantics(grid_base3, grididx_base3) == expected_quantics_base3
    @test quantics_to_grididx(grid_base3, expected_quantics_base3) == grididx_base3
end

@testitem "quantics_to_grididx ∘ grididx_to_quantics == identity" begin
    grid = DiscretizedGrid((5, 3, 17); base=13)
    for _ in 1:100
        grididx = ntuple(d -> rand(1:QuanticsGrids.grid_Rs(grid)[d]), ndims(grid))
        @test quantics_to_grididx(grid, grididx_to_quantics(grid, grididx)) == grididx
    end
end

@testitem "grididx_to_quantics ∘ quantics_to_grididx == identity" begin
    grid = DiscretizedGrid((48, 31, 62))
    for _ in 1:100
        quantics = rand(1:2, length(grid))
        @test grididx_to_quantics(grid, quantics_to_grididx(grid, quantics)) == quantics
    end
end

@testitem "grididx_to_quantics ∘ quantics_to_grididx == identity, base != 2" begin
    base = 7
    grid = DiscretizedGrid((22, 9, 14); base)
    for _ in 1:100
        quantics = rand(1:base, length(grid))
        @test grididx_to_quantics(grid, quantics_to_grididx(grid, quantics)) == quantics
    end
end

@testitem "ctor from indextable" begin
    grid = DiscretizedGrid((:a, :b, :c), [[(:a, 4)], [(:a, 3)], [(:a, 2)], [(:a, 1)], [(:b, 1)], [(:b, 2)], [(:b, 3)], [(:c, 1)], [(:c, 2)], [(:c, 3)]])
    @test QuanticsGrids.grid_Rs(grid) == (4, 3, 3)
    @test QuanticsGrids.lower_bound(grid) == (0.0, 0.0, 0.0)
    @test QuanticsGrids.upper_bound(grid) == (1.0, 1.0, 1.0)
    @test grid.discretegrid.variablenames == (:a, :b, :c)
end

@testitem "ctor from indextable, quantics <-> grididx" begin
    grid = DiscretizedGrid((:a, :b, :c), [[(:a, 4)], [(:a, 3)], [(:a, 2)], [(:a, 1)], [(:b, 1)], [(:b, 2)], [(:b, 3)], [(:c, 1)], [(:c, 2)], [(:c, 3)]])
    @test quantics_to_grididx(grid, [1, 2, 1, 2, 1, 2, 1, 2, 1, 2]) == (11, 3, 6)
    @test grididx_to_quantics(grid, (11, 3, 6)) == [1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
end

@testitem "ctor from indextable, quantics <-> grididx, fused indices" begin
    grid = DiscretizedGrid((:a, :b, :c, :d), [[(:a, 4)], [(:a, 3)], [(:a, 2)], [(:a, 1)], [(:b, 1), (:d, 1)], [(:b, 2)], [(:b, 3)], [(:c, 1), (:d, 2)], [(:c, 2)], [(:c, 3)]])
    @test quantics_to_grididx(grid, [1, 2, 1, 2, 2, 2, 1, 4, 1, 2]) == (11, 3, 6, 4)
    @test grididx_to_quantics(grid, (11, 3, 6, 4)) == [1, 2, 1, 2, 2, 2, 1, 4, 1, 2]
end

@testitem "challenging tests - extreme edge cases" begin
    # Test with minimum valid grididx (all 1s)
    grid = DiscretizedGrid((10, 5, 8); unfoldingscheme=:interleaved)
    min_grididx = (1, 1, 1)
    quantics = grididx_to_quantics(grid, min_grididx)
    @test all(q -> q == 1, quantics)
    @test quantics_to_grididx(grid, quantics) == min_grididx

    # Test with maximum valid grididx
    max_grididx = (2^10, 2^5, 2^8)
    quantics = grididx_to_quantics(grid, max_grididx)
    @test all(q -> q == 2, quantics)
    @test quantics_to_grididx(grid, quantics) == max_grididx
end

@testitem "challenging tests - mixed bases" begin
    # Test with base 3
    grid = DiscretizedGrid((4, 6, 3); base=3)
    for _ in 1:50
        quantics = rand(1:3, length(grid))
        grididx = quantics_to_grididx(grid, quantics)
        @test grididx_to_quantics(grid, grididx) == quantics
        @test all(1 .<= grididx .<= (3 .^ QuanticsGrids.grid_Rs(grid)))
    end

    # Test with base 5
    grid = DiscretizedGrid((3, 4); base=5, unfoldingscheme=:interleaved)
    for _ in 1:50
        grididx = (rand(1:5^3), rand(1:5^4))
        quantics = grididx_to_quantics(grid, grididx)
        @test quantics_to_grididx(grid, quantics) == grididx
        @test all(1 .<= quantics .<= 5)
    end
end

@testitem "challenging tests - complex fused indices" begin
    # Multiple variables fused in single sites
    grid = DiscretizedGrid(
        (:x, :y, :z),
        [
            [(:x, 3), (:y, 2), (:z, 1)],  # 3 variables in one site
            [(:x, 2)],                     # single variable
            [(:y, 1), (:z, 2)],          # 2 variables fused
            [(:x, 1)],                     # single variable
            [(:z, 3)]                      # single variable
        ]
    )

    # Test specific known values
    @test quantics_to_grididx(grid, [8, 1, 4, 2, 1]) == (6, 4, 7)
    @test grididx_to_quantics(grid, (6, 4, 7)) == [8, 1, 4, 2, 1]

    # Test random values
    for _ in 1:100
        quantics = [
            rand(1:8),   # site 1: base^3 = 8 possibilities
            rand(1:2),   # site 2: base = 2 possibilities  
            rand(1:4),   # site 3: base^2 = 4 possibilities
            rand(1:2),   # site 4: base = 2 possibilities
            rand(1:2)    # site 5: base = 2 possibilities
        ]
        grididx = quantics_to_grididx(grid, quantics)
        @test grididx_to_quantics(grid, grididx) == quantics
    end
end

@testitem "challenging tests - asymmetric grids" begin
    # Very asymmetric grid dimensions
    grid = DiscretizedGrid((20, 0, 15, 3))

    for _ in 1:50
        # Test boundary conditions
        quantics = rand(1:2, length(grid))
        grididx = quantics_to_grididx(grid, quantics)
        @test grididx_to_quantics(grid, grididx) == quantics

        # Verify grid index bounds
        @test 1 <= grididx[1] <= 2^20
        @test grididx[2] == 1  # R=0 means only one possible value
        @test 1 <= grididx[3] <= 2^15
        @test 1 <= grididx[4] <= 2^3
    end
end

@testitem "challenging tests - single dimension edge cases" begin
    # 1D grid with large R
    grid = DiscretizedGrid((25,))

    # Test extremes
    min_quantics = ones(Int, 25)
    max_quantics = fill(2, 25)

    @test quantics_to_grididx(grid, min_quantics) == 1
    @test quantics_to_grididx(grid, max_quantics) == 2^25
    @test grididx_to_quantics(grid, (1,)) == min_quantics
    @test grididx_to_quantics(grid, (2^25,)) == max_quantics

    # Test middle values
    for _ in 1:20
        quantics = rand(1:2, 25)
        grididx = quantics_to_grididx(grid, quantics)
        @test grididx_to_quantics(grid, grididx) == quantics
    end
end

@testitem "challenging tests - high dimensional grids" begin
    # 8D grid with moderate R values
    grid = DiscretizedGrid(ntuple(i -> 4 + (i % 3), 8); base=3)

    for _ in 1:30
        quantics = rand(1:3, length(grid))
        grididx = quantics_to_grididx(grid, quantics)
        @test grididx_to_quantics(grid, grididx) == quantics

        # Verify all grid indices are within bounds
        for d in 1:8
            @test 1 <= grididx[d] <= 3^QuanticsGrids.grid_Rs(grid)[d]
        end
    end
end

@testitem "challenging tests - stress test with complex patterns" begin
    # Grid with intentionally complex index table structure
    grid = DiscretizedGrid(
        (:a, :b, :c, :d, :e),
        [
            [(:e, 1)],                           # site 1
            [(:a, 5), (:c, 4)],                 # site 2: 2 vars fused
            [(:b, 3)],                           # site 3
            [(:a, 4), (:b, 2), (:d, 3)],       # site 4: 3 vars fused
            [(:c, 3), (:e, 2)],                 # site 5: 2 vars fused
            [(:a, 3)],                           # site 6
            [(:b, 1), (:d, 2)],                 # site 7: 2 vars fused  
            [(:a, 2), (:c, 2), (:d, 1), (:e, 3)], # site 8: 4 vars fused
            [(:a, 1)],                           # site 9
            [(:c, 1)]                            # site 10
        ];
        base=3
    )

    # Stress test with many random conversions
    for _ in 1:200
        # Generate valid quantics vector
        quantics = [
            rand(1:3),    # site 1
            rand(1:9),    # site 2: 3^2 = 9
            rand(1:3),    # site 3
            rand(1:27),   # site 4: 3^3 = 27
            rand(1:9),    # site 5: 3^2 = 9
            rand(1:3),    # site 6
            rand(1:9),    # site 7: 3^2 = 9
            rand(1:81),   # site 8: 3^4 = 81
            rand(1:3),    # site 9
            rand(1:3)     # site 10
        ]

        grididx = quantics_to_grididx(grid, quantics)
        recovered = grididx_to_quantics(grid, grididx)
        @test recovered == quantics

        # Verify grid indices are reasonable
        @test all(1 .<= grididx .<= (3 .^ QuanticsGrids.grid_Rs(grid)))
    end
end

@testitem "fused ordering compatibility with DiscretizedGrid" begin
    using QuanticsGrids: DiscretizedGrid

    # Test 2D case - this was the main issue discovered during benchmarking
    # where coordinates were being swapped between implementations
    R = 3  # Use same R for both dimensions to create equivalent grids
    lower = (0.0, 0.0)
    upper = (1.0, 1.0)

    # Create equivalent grids
    old_grid = DiscretizedGrid{2}(R, lower, upper; unfoldingscheme=:fused)
    new_grid = DiscretizedGrid((R, R); lower_bound=lower, upper_bound=upper, unfoldingscheme=:fused)

    # Test specific quantics vectors that revealed the ordering issue
    test_cases = [
        [1, 2, 1],  # This was giving different coordinates before the fix
        [2, 1, 2],
        [1, 1, 1],  # Edge case: all 1s
        [2, 2, 2],  # Edge case: all 2s
        [2, 1, 1],
        [1, 2, 2]
    ]

    for quantics in test_cases
        # Convert to coordinates using both implementations
        old_coord = quantics_to_origcoord(old_grid, quantics)
        new_coord = quantics_to_origcoord(new_grid, quantics)

        # Coordinates should be identical (within floating point precision)
        @test all(old_coord .≈ new_coord)

        # Also verify the round-trip consistency within each implementation
        old_grididx = quantics_to_grididx(old_grid, quantics)
        new_grididx = quantics_to_grididx(new_grid, quantics)
        @test old_grididx == new_grididx

        # Verify quantics -> grididx -> origcoord is consistent
        old_coord_via_grididx = grididx_to_origcoord(old_grid, old_grididx)
        new_coord_via_grididx = grididx_to_origcoord(new_grid, new_grididx)
        @test all(old_coord_via_grididx .≈ new_coord_via_grididx)
    end

    # Test 3D case to ensure the fix works for higher dimensions
    R_3d = 2  # Use same R for all dimensions
    lower_3d = (0.0, 0.0, 0.0)
    upper_3d = (1.0, 1.0, 1.0)

    old_grid_3d = DiscretizedGrid{3}(R_3d, lower_3d, upper_3d; unfoldingscheme=:fused)
    new_grid_3d = DiscretizedGrid((R_3d, R_3d, R_3d); lower_bound=lower_3d, upper_bound=upper_3d, unfoldingscheme=:fused)

    # Test several 3D quantics vectors
    test_cases_3d = [
        [1, 1],  # Length R=2 for fused scheme
        [2, 2],
        [1, 2],
        [2, 1]
    ]

    for quantics in test_cases_3d
        old_coord = quantics_to_origcoord(old_grid_3d, quantics)
        new_coord = quantics_to_origcoord(new_grid_3d, quantics)
        @test all(old_coord .≈ new_coord)
    end
end
