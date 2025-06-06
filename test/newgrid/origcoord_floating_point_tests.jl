@testitem "floating point precision - coordinates exactly at boundaries" begin
    using QuanticsGrids
    # Test behavior when coordinates are exactly at grid boundaries
    R = 10
    grid = NewDiscretizedGrid{1}(R)

    # Test exact lower bound
    @test QuanticsGrids.origcoord_to_grididx(grid, 0.0) == 1

    # Test exact upper bound (should be at grid_max, not upper_bound)
    grid_max_val = QuanticsGrids.grid_max(grid)
    @test QuanticsGrids.origcoord_to_grididx(grid, grid_max_val) == 2^R

    # Test coordinates that are exactly representable as grid points
    step = QuanticsGrids.grid_step(grid)
    for i in 1:min(100, 2^R)
        exact_coord = (i - 1) * step
        @test QuanticsGrids.origcoord_to_grididx(grid, exact_coord) == i
    end
end

@testitem "floating point precision - very small step sizes" begin
    using QuanticsGrids
    # Test with large R values that create very small step sizes
    R = 50  # Creates step size of 1/2^50 ≈ 8.9e-16
    grid = NewDiscretizedGrid{1}(R)

    step = QuanticsGrids.grid_step(grid)
    @test step ≈ 1.0 / 2^R

    # Test coordinates that should map to specific grid indices
    # Use Float64 arithmetic carefully to avoid precision loss
    test_indices = [1, 2, 3, 2^10, 2^20, 2^R ÷ 2, 2^R - 1, 2^R]

    for idx in test_indices
        # Calculate coordinate using exact arithmetic when possible
        coord = (idx - 1) * step
        recovered_idx = QuanticsGrids.origcoord_to_grididx(grid, coord)
        @test recovered_idx == idx
    end
end

@testitem "floating point precision - catastrophic cancellation" begin
    using QuanticsGrids
    # Test cases where subtraction coordinate - lower_bound might lose precision

    # Case 1: Very close lower_bound and coordinate values
    lower = 1e15
    upper = 1e15 + 1.0
    grid = NewDiscretizedGrid{1}(10; lower_bound=(lower,), upper_bound=(upper,))

    # Test coordinates very close to lower bound
    eps_values = [1e-15, 1e-14, 1e-13, 1e-12, 1e-10]
    for eps in eps_values
        coord = lower + eps
        if coord <= QuanticsGrids.upper_bound(grid)
            idx = QuanticsGrids.origcoord_to_grididx(grid, coord)
            @test idx >= 1 && idx <= 2^10

            # Test round-trip
            recovered_coord = grididx_to_origcoord(grid, idx)
            recovered_idx = QuanticsGrids.origcoord_to_grididx(grid, recovered_coord)
            @test recovered_idx == idx
        end
    end

    # Case 2: Very large coordinate values where precision is lost
    lower = 1e12
    upper = 1e12 + 1e6
    grid2 = NewDiscretizedGrid{1}(15; lower_bound=(lower,), upper_bound=(upper,))

    step = QuanticsGrids.grid_step(grid2)
    for i in [1, 100, 1000, 2^15 ÷ 2, 2^15]
        coord = lower + (i - 1) * step
        if coord <= QuanticsGrids.upper_bound(grid2)
            idx = QuanticsGrids.origcoord_to_grididx(grid2, coord)
            @test abs(idx - i) <= 1
        end
    end
end

@testitem "floating point precision - round-trip consistency stress test" begin
    using QuanticsGrids
    # Test round-trip consistency: grididx -> origcoord -> grididx
    test_grids = [
        # Standard precision range
        NewDiscretizedGrid{1}(20),
        # Very small range
        NewDiscretizedGrid{1}(15; lower_bound=(1e-10,), upper_bound=(1e-9,)),
        # Very large range  
        NewDiscretizedGrid{1}(12; lower_bound=(1e10,), upper_bound=(1e11,)),
        # Negative range
        NewDiscretizedGrid{1}(18; lower_bound=(-1e5,), upper_bound=(-1e4,)),
        # Asymmetric around zero
        NewDiscretizedGrid{1}(16; lower_bound=(-1e-6,), upper_bound=(1e-5,)),
    ]

    for grid in test_grids
        max_idx = 2^grid.Rs[1]
        # Test a representative sample of indices including edge cases
        test_indices = unique([1, 2, 3, max_idx ÷ 4, max_idx ÷ 2, 3 * max_idx ÷ 4, max_idx - 2, max_idx - 1, max_idx])

        for idx in test_indices
            # Forward: grididx -> origcoord -> grididx
            coord = grididx_to_origcoord(grid, idx)
            recovered_idx = QuanticsGrids.origcoord_to_grididx(grid, coord)
            @test recovered_idx == idx
        end
    end
end

@testitem "floating point precision - coordinates near but not exactly at boundaries" begin
    using QuanticsGrids
    R = 15
    grid = NewDiscretizedGrid{1}(R)

    # Test coordinates that are very close to boundaries but not exactly at them
    eps_values = [eps(Float64), 2 * eps(Float64), 1e-15, 1e-14, 1e-13, 1e-12]

    lower = QuanticsGrids.lower_bound(grid)
    upper = QuanticsGrids.upper_bound(grid)

    for eps_val in eps_values
        # Just above lower bound
        coord_above_lower = lower + eps_val
        if coord_above_lower <= upper
            idx = QuanticsGrids.origcoord_to_grididx(grid, coord_above_lower)
            @test idx >= 1
        end

        # Just below upper bound
        coord_below_upper = upper - eps_val
        if coord_below_upper >= lower
            idx = QuanticsGrids.origcoord_to_grididx(grid, coord_below_upper)
            @test idx <= 2^R
        end
    end
end

@testitem "floating point precision - multidimensional edge cases" begin
    using QuanticsGrids
    # Test 2D grids where floating point errors might compound
    R = 12
    grid = NewDiscretizedGrid((R, R);
        lower_bound=(0.0, -1.0),
        upper_bound=(1.0, 1.0))

    step = QuanticsGrids.grid_step(grid)
    max_indices = (2^R, 2^R)

    # Test corner cases
    corner_indices = [(1, 1), (1, max_indices[2]), (max_indices[1], 1), max_indices]

    for corner_idx in corner_indices
        coord = grididx_to_origcoord(grid, corner_idx)
        recovered_idx = QuanticsGrids.origcoord_to_grididx(grid, coord)
        @test recovered_idx == corner_idx
    end

    # Test coordinates computed with potential precision issues
    for _ in 1:50
        # Generate random indices
        idx = (rand(1:2^R), rand(1:2^R))

        # Convert to coordinates and back
        coord = grididx_to_origcoord(grid, idx)
        recovered_idx = QuanticsGrids.origcoord_to_grididx(grid, coord)

        @test recovered_idx == idx
    end
end

@testitem "floating point precision - extreme coordinate ranges" begin
    using QuanticsGrids
    # Test grids with extreme coordinate ranges that stress floating point arithmetic

    test_cases = [
        # Very small positive range
        (1e-16, 1e-15, 10),
        # Very large range
        (1e15, 1e16, 8),
        # Range crossing zero with very small values
        (-1e-15, 1e-15, 12),
        # Asymmetric range with vastly different magnitudes
        (-1e10, 1e-8, 10),
        # Range where step size is close to machine epsilon
        (0.0, 1e-14, 8),
    ]

    for (lower, upper, R) in test_cases
        grid = NewDiscretizedGrid{1}(R; lower_bound=(lower,), upper_bound=(upper,))

        # Test boundary indices
        boundary_indices = [1, 2, 2^R - 1, 2^R]

        for idx in boundary_indices
            coord = grididx_to_origcoord(grid, idx)

            # Check that coordinate is within expected bounds
            @test QuanticsGrids.lower_bound(grid) <= coord <= QuanticsGrids.upper_bound(grid)

            # Test round-trip
            recovered_idx = QuanticsGrids.origcoord_to_grididx(grid, coord)
            @test recovered_idx == idx
        end
    end
end

@testitem "floating point precision - step size edge cases" begin
    using QuanticsGrids
    # Test cases where step size computation might be problematic

    # Case 1: Step size that's not exactly representable in Float64
    grid1 = NewDiscretizedGrid{1}(10; lower_bound=(0.0,), upper_bound=(1.0 / 3.0,))
    step1 = QuanticsGrids.grid_step(grid1)

    # Test that step size calculations are consistent
    @test step1 ≈ (1.0 / 3.0) / 2^10

    for i in 1:2^10
        coord = grididx_to_origcoord(grid1, i)
        recovered_idx = QuanticsGrids.origcoord_to_grididx(grid1, coord)
        @test recovered_idx == i
    end

    # Case 2: Step size involving irrational numbers
    sqrt2 = sqrt(2.0)
    grid2 = NewDiscretizedGrid{1}(8; lower_bound=(0.0,), upper_bound=(sqrt2,))

    for i in [1, 2^4, 2^8]
        coord = grididx_to_origcoord(grid2, i)
        recovered_idx = QuanticsGrids.origcoord_to_grididx(grid2, coord)
        @test recovered_idx == i
    end

    # Case 3: Very small step size that approaches machine epsilon
    grid3 = NewDiscretizedGrid{1}(4; lower_bound=(0.0,), upper_bound=(eps(Float64) * 100,))

    for i in 1:2^4
        coord = grididx_to_origcoord(grid3, i)
        recovered_idx = QuanticsGrids.origcoord_to_grididx(grid3, coord)
        @test recovered_idx == i
    end
end

@testitem "floating point precision - boundary coordinate generation" begin
    using QuanticsGrids
    # Test coordinates that are generated by arithmetic operations and might not be exact
    R = 20
    grid = NewDiscretizedGrid{1}(R)

    step = QuanticsGrids.grid_step(grid)

    # Test coordinates generated by accumulative addition (potential error accumulation)
    coord = 0.0
    for i in 1:min(1000, 2^R)
        expected_idx = i
        actual_idx = QuanticsGrids.origcoord_to_grididx(grid, coord)

        # Allow for small errors due to accumulation
        @test abs(actual_idx - expected_idx) <= 1

        global coord += step
        if coord > QuanticsGrids.upper_bound(grid)
            break
        end
    end
end

@testitem "floating point precision - clamp behavior verification" begin
    using QuanticsGrids
    # Test that the clamping behavior works correctly for edge cases
    R = 8
    grid = NewDiscretizedGrid{1}(R)

    # Test coordinates that would round to indices outside valid range
    step = QuanticsGrids.grid_step(grid)

    # Coordinate slightly before first grid point (should clamp to 1)
    coord_before = -step / 2
    if coord_before >= QuanticsGrids.lower_bound(grid)
        idx = QuanticsGrids.origcoord_to_grididx(grid, coord_before)
        @test idx == 1
    end

    # Coordinate slightly after last grid point (should clamp to max index)
    last_coord = grididx_to_origcoord(grid, 2^R)
    coord_after = last_coord + step / 2
    if coord_after <= QuanticsGrids.upper_bound(grid)
        idx = QuanticsGrids.origcoord_to_grididx(grid, coord_after)
        @test idx == 2^R
    end

    # Test that coordinates exactly at grid_max map to the last index
    grid_max_coord = QuanticsGrids.grid_max(grid)
    idx_at_max = QuanticsGrids.origcoord_to_grididx(grid, grid_max_coord)
    @test idx_at_max == 2^R
end

@testitem "floating point precision - includeendpoint behavior edge cases" begin
    using QuanticsGrids
    # Test floating point behavior with includeendpoint=true vs false
    R = 10

    grid_without = NewDiscretizedGrid{1}(R; includeendpoint=false)
    grid_with = NewDiscretizedGrid{1}(R; includeendpoint=true)

    # Test that the upper bound is handled differently
    upper_bound = 1.0

    # Without includeendpoint: coordinate at upper_bound should be valid but map to last index
    if upper_bound <= QuanticsGrids.upper_bound(grid_without)
        idx_without = QuanticsGrids.origcoord_to_grididx(grid_without, upper_bound)
        @test idx_without == 2^R
    end

    # With includeendpoint: coordinate at upper_bound should map to last index
    idx_with = QuanticsGrids.origcoord_to_grididx(grid_with, upper_bound)
    @test idx_with == 2^R

    # Test round-trip consistency for both cases
    for grid in [grid_without, grid_with]
        step = QuanticsGrids.grid_step(grid)

        # Test coordinates near the upper boundary
        test_coords = [
            QuanticsGrids.grid_max(grid),
            QuanticsGrids.grid_max(grid) - step / 10,
            QuanticsGrids.grid_max(grid) + step / 10,  # This might be out of bounds
        ]

        for coord in test_coords
            if QuanticsGrids.lower_bound(grid) <= coord <= QuanticsGrids.upper_bound(grid)
                idx = QuanticsGrids.origcoord_to_grididx(grid, coord)
                recovered_coord = grididx_to_origcoord(grid, idx)
                recovered_idx = QuanticsGrids.origcoord_to_grididx(grid, recovered_coord)
                @test recovered_idx == idx
            end
        end
    end
end

@testitem "floating point precision - mixed precision operations" begin
    using QuanticsGrids
    # Test scenarios where different precision operations might interact

    # Create grid with parameters that don't divide evenly
    grid = NewDiscretizedGrid{1}(12; lower_bound=(0.1,), upper_bound=(0.9,))

    step = QuanticsGrids.grid_step(grid)

    # Test coordinates calculated using different methods
    for i in [1, 100, 1000, 2^12]
        # Method 1: Direct calculation
        coord1 = 0.1 + (i - 1) * step

        # Method 2: Via grid function
        coord2 = grididx_to_origcoord(grid, i)

        # Both should give same results (within floating point precision)
        @test coord1 ≈ coord2

        # Both should map back to same index
        if coord1 <= QuanticsGrids.upper_bound(grid)
            idx1 = QuanticsGrids.origcoord_to_grididx(grid, coord1)
            @test idx1 == i
        end

        idx2 = QuanticsGrids.origcoord_to_grididx(grid, coord2)
        @test idx2 == i
    end
end
