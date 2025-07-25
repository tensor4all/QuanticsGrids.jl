@testitem "InherentDiscreteGrid" begin
    for unfoldingscheme in [:interleaved, :fused],
        step in [(1, 1, 1), (1, 1, 2)],
        origin in [(1, 1, 1), (1, 1, 2)]

        m = InherentDiscreteGrid{3}(5, origin; unfoldingscheme, step=step)
        @test QuanticsGrids.grid_min(m) == origin
        @test QuanticsGrids.grid_step(m) == step

        if unfoldingscheme === :interleaved
            @test QuanticsGrids.localdimensions(m) == fill(2, 3 * 5)
        else
            @test QuanticsGrids.localdimensions(m) == fill(2^3, 5)
        end

        for idx in [(1, 1, 1), (1, 1, 2), (1, 25, 1), (14, 1, 1), (25, 25, 25)]
            c = QuanticsGrids.grididx_to_origcoord(m, idx)
            @test QuanticsGrids.origcoord_to_grididx(m, c) == idx

            q = QuanticsGrids.grididx_to_quantics(m, idx)
            if unfoldingscheme === :fused
                @test length(q) == 5
            else
                @test length(q) == 3 * 5
            end
            @test all((1 .<= q) .&& (q .<= 2^3))
            @test QuanticsGrids.quantics_to_origcoord(m, q) == c
        end
    end
end

@testitem "InherentDiscreteGrid constructors" begin
    # Test basic 1D constructor
    R = 4
    grid_1d = InherentDiscreteGrid{1}(R, 1; base=2, step=1, unfoldingscheme=:fused)
    @test grid_1d.base == 2
    @test grid_1d.origin == (1,)
    @test grid_1d.step == (1,)

    # Test multi-dimensional constructor with tuple origin and step
    origin_3d = (5, 10, 15)
    step_3d = (2, 3, 1)
    grid_3d = InherentDiscreteGrid{3}(R, origin_3d; base=3, step=step_3d, unfoldingscheme=:interleaved)
    @test grid_3d.base == 3
    @test grid_3d.origin == origin_3d
    @test grid_3d.step == step_3d

    # Test default parameters
    grid_default = InherentDiscreteGrid{2}(R, (0, 0))
    @test grid_default.base == 2
    @test grid_default.step == (1, 1)
end

@testitem "InherentDiscreteGrid convenience constructors" begin
    # Test single integer origin constructor
    R = 3
    grid_single = InherentDiscreteGrid{2}(R, 5; step=2)
    @test grid_single.origin == (5, 5)
    @test grid_single.step == (2, 2)

    # Test mixed single/tuple parameters
    grid_mixed = InherentDiscreteGrid{3}(R, (1, 2, 3); step=4)
    @test grid_mixed.origin == (1, 2, 3)
    @test grid_mixed.step == (4, 4, 4)
end

@testitem "InherentDiscreteGrid basic properties" begin
    R = 5
    base = 3
    origin = (10, 20)
    step = (2, 5)

    grid = InherentDiscreteGrid{2}(R, origin; base=base, step=step)

    # Test basic properties
    @test ndims(grid) == 2
    @test length(grid) == R  # Number of quantics sites for fused
    @test grid.base == base

    # Test grid boundaries
    @test QuanticsGrids.grid_min(grid) == origin
    @test QuanticsGrids.grid_max(grid) == origin .+ step .* (base^R - 1)
    @test QuanticsGrids.grid_step(grid) == step
    @test QuanticsGrids.grid_origin(grid) == origin

    # Test with interleaved scheme
    grid_interleaved = InherentDiscreteGrid{2}(R, origin; base=base, step=step, unfoldingscheme=:interleaved)
    @test length(grid_interleaved) == 2 * R  # Number of quantics sites for interleaved
end

@testitem "InherentDiscreteGrid coordinate conversions" begin
    R = 4
    base = 2
    origin = (5, 10, 15)
    step = (2, 3, 1)

    grid = InherentDiscreteGrid{3}(R, origin; base=base, step=step)

    # Test grididx_to_origcoord and origcoord_to_grididx roundtrip
    test_grididx = [(1, 1, 1), (1, 1, 2), (1, 8, 1), (4, 1, 1), (16, 16, 16)]
    for grididx in test_grididx
        origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
        @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx
    end

    # Test origcoord calculation formula
    grididx = (3, 5, 7)
    expected_origcoord = origin .+ step .* (grididx .- 1)
    @test QuanticsGrids.grididx_to_origcoord(grid, grididx) == expected_origcoord

    # Test boundary cases
    min_grididx = (1, 1, 1)
    max_grididx = (base^R, base^R, base^R)
    @test QuanticsGrids.grididx_to_origcoord(grid, min_grididx) == origin
    @test QuanticsGrids.grididx_to_origcoord(grid, max_grididx) == origin .+ step .* (base^R - 1)
end

@testitem "InherentDiscreteGrid quantics conversions" begin
    # Test fused quantics
    R = 3
    base = 2
    grid_fused = InherentDiscreteGrid{2}(R, (1, 1); base=base, unfoldingscheme=:fused)

    # Test grididx_to_quantics and quantics_to_grididx roundtrip
    test_grididx = [(1, 1), (1, 2), (2, 1), (4, 8), (8, 8)]
    for grididx in test_grididx
        quantics = QuanticsGrids.grididx_to_quantics(grid_fused, grididx)
        @test QuanticsGrids.quantics_to_grididx(grid_fused, quantics) == grididx
        @test length(quantics) == R  # Fused should have R quantics
        @test all(1 .<= quantics .<= base^2)  # Each quantics should be in valid range
    end

    # Test interleaved quantics
    grid_interleaved = InherentDiscreteGrid{2}(R, (1, 1); base=base, unfoldingscheme=:interleaved)

    for grididx in test_grididx
        quantics = QuanticsGrids.grididx_to_quantics(grid_interleaved, grididx)
        @test QuanticsGrids.quantics_to_grididx(grid_interleaved, quantics) == grididx
        @test length(quantics) == 2 * R  # Interleaved should have 2*R quantics
        @test all(1 .<= quantics .<= base)  # Each quantics should be in valid range
    end
end

@testitem "InherentDiscreteGrid quantics_to_origcoord and origcoord_to_quantics" begin
    R = 4
    origin = (3, 7)
    step = (2, 3)
    grid = InherentDiscreteGrid{2}(R, origin; step=step)

    # Test quantics_to_origcoord and origcoord_to_quantics roundtrip
    for _ in 1:20
        # Generate random valid quantics
        quantics = rand(1:2^2, R)  # Each quantics element should be valid for the base and dimension

        origcoord = QuanticsGrids.quantics_to_origcoord(grid, quantics)
        @test QuanticsGrids.origcoord_to_quantics(grid, origcoord) == quantics

        # Verify the conversion chain: quantics -> grididx -> origcoord
        grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
        expected_origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
        @test origcoord == expected_origcoord
    end
end

@testitem "InherentDiscreteGrid localdimensions" begin
    R = 4
    base = 3

    # Test fused scheme
    grid_fused_1d = InherentDiscreteGrid{1}(R, 1; base=base, unfoldingscheme=:fused)
    @test QuanticsGrids.localdimensions(grid_fused_1d) == fill(base, R)

    grid_fused_3d = InherentDiscreteGrid{3}(R, (1, 1, 1); base=base, unfoldingscheme=:fused)
    @test QuanticsGrids.localdimensions(grid_fused_3d) == fill(base^3, R)

    # Test interleaved scheme
    grid_interleaved_1d = InherentDiscreteGrid{1}(R, 1; base=base, unfoldingscheme=:interleaved)
    @test QuanticsGrids.localdimensions(grid_interleaved_1d) == fill(base, R)

    grid_interleaved_3d = InherentDiscreteGrid{3}(R, (1, 1, 1); base=base, unfoldingscheme=:interleaved)
    @test QuanticsGrids.localdimensions(grid_interleaved_3d) == fill(base, 3 * R)

    # Test different bases
    for test_base in [2, 4, 5]
        grid = InherentDiscreteGrid{2}(R, (1, 1); base=test_base, unfoldingscheme=:fused)
        @test QuanticsGrids.localdimensions(grid) == fill(test_base^2, R)
    end
end

@testitem "InherentDiscreteGrid different bases" begin
    R = 3
    origin = (1, 1)

    for base in [2, 3, 4, 5]
        grid = InherentDiscreteGrid{2}(R, origin; base=base)

        # Test that grid can handle different bases
        max_grididx = (base^R, base^R)
        min_grididx = (1, 1)

        # Test coordinate conversion works with different bases
        max_origcoord = QuanticsGrids.grididx_to_origcoord(grid, max_grididx)
        min_origcoord = QuanticsGrids.grididx_to_origcoord(grid, min_grididx)

        @test QuanticsGrids.origcoord_to_grididx(grid, max_origcoord) == max_grididx
        @test QuanticsGrids.origcoord_to_grididx(grid, min_origcoord) == min_grididx

        # Test quantics conversion with different bases
        quantics = QuanticsGrids.grididx_to_quantics(grid, (2, 3))
        @test all(1 .<= quantics .<= base^2)  # For fused scheme
        @test QuanticsGrids.quantics_to_grididx(grid, quantics) == (2, 3)
    end
end

@testitem "InherentDiscreteGrid scalar vs tuple return values" begin
    R = 3

    # Test 1D grid returns scalars
    grid_1d = InherentDiscreteGrid{1}(R, 5; step=2)
    grididx_1d = 3
    origcoord_1d = QuanticsGrids.grididx_to_origcoord(grid_1d, grididx_1d)
    @test origcoord_1d isa Int
    @test QuanticsGrids.origcoord_to_grididx(grid_1d, origcoord_1d) == grididx_1d

    @test QuanticsGrids.grid_min(grid_1d) isa Int
    @test QuanticsGrids.grid_max(grid_1d) isa Int
    @test QuanticsGrids.grid_step(grid_1d) isa Int
    @test QuanticsGrids.grid_origin(grid_1d) isa Int

    # Test multi-dimensional grid returns tuples
    grid_2d = InherentDiscreteGrid{2}(R, (3, 7); step=(1, 2))
    grididx_2d = (2, 4)
    origcoord_2d = QuanticsGrids.grididx_to_origcoord(grid_2d, grididx_2d)
    @test origcoord_2d isa Tuple{Int,Int}
    @test QuanticsGrids.origcoord_to_grididx(grid_2d, origcoord_2d) == grididx_2d

    @test QuanticsGrids.grid_min(grid_2d) isa Tuple{Int,Int}
    @test QuanticsGrids.grid_max(grid_2d) isa Tuple{Int,Int}
    @test QuanticsGrids.grid_step(grid_2d) isa Tuple{Int,Int}
    @test QuanticsGrids.grid_origin(grid_2d) isa Tuple{Int,Int}
end

@testitem "InherentDiscreteGrid comprehensive conversion test" begin
    # Test all possible conversion combinations like in the old grid tests
    R = 4
    base = 3
    origin = (2, 5)
    step = (3, 2)

    grid = InherentDiscreteGrid{2}(R, origin; base=base, step=step)

    reprs = [:grididx, :quantics, :origcoord]

    transforms = Dict(
        (:grididx, :quantics) => QuanticsGrids.grididx_to_quantics,
        (:grididx, :origcoord) => QuanticsGrids.grididx_to_origcoord,
        (:quantics, :grididx) => QuanticsGrids.quantics_to_grididx,
        (:quantics, :origcoord) => QuanticsGrids.quantics_to_origcoord,
        (:origcoord, :grididx) => QuanticsGrids.origcoord_to_grididx,
        (:origcoord, :quantics) => QuanticsGrids.origcoord_to_quantics,
    )

    initial_grididx = (2, 3)
    data = Dict{Symbol,Any}()
    data[:grididx] = initial_grididx

    # Build up all representations
    while true
        newdata = false
        for src in reprs, dst in reprs
            if src == dst
                continue
            end
            if !(src ∈ keys(data))
                continue
            end

            if dst ∈ keys(data)
                @test transforms[(src, dst)](grid, data[src]) == data[dst]
            else
                data[dst] = transforms[(src, dst)](grid, data[src])
                newdata = true
            end
        end
        if newdata == false
            break
        end
    end

    # Verify all conversions are consistent
    @test length(data) == 3
    @test haskey(data, :grididx)
    @test haskey(data, :quantics)
    @test haskey(data, :origcoord)
end

@testitem "InherentDiscreteGrid edge cases and boundary conditions" begin
    R = 3
    base = 2

    # Test minimum grid size (R=1)
    grid_min = InherentDiscreteGrid{1}(1, 10; base=base, step=5)
    @test QuanticsGrids.grididx_to_origcoord(grid_min, 1) == 10
    @test QuanticsGrids.grididx_to_origcoord(grid_min, base) == 10 + 5 * (base - 1)

    # Test large step sizes
    grid_large_step = InherentDiscreteGrid{2}(R, (0, 0); step=(100, 200))
    grididx = (2, 3)
    expected_coord = (100, 400)  # (0 + 100*(2-1), 0 + 200*(3-1))
    @test QuanticsGrids.grididx_to_origcoord(grid_large_step, grididx) == expected_coord

    # Test negative origins
    grid_neg = InherentDiscreteGrid{2}(R, (-10, -20); step=(2, 3))
    min_coord = QuanticsGrids.grididx_to_origcoord(grid_neg, (1, 1))
    @test min_coord == (-10, -20)

    # Test single dimension
    grid_1d = InherentDiscreteGrid{1}(R, 5; step=3)
    for i in 1:base^R
        coord = QuanticsGrids.grididx_to_origcoord(grid_1d, i)
        @test QuanticsGrids.origcoord_to_grididx(grid_1d, coord) == i
    end
end

@testitem "InherentDiscreteGrid consistency with unfolding schemes" begin
    R = 4
    origin = (1, 1, 1)
    step = (1, 2, 3)
    base = 2

    grid_fused = InherentDiscreteGrid{3}(R, origin; step=step, base=base, unfoldingscheme=:fused)
    grid_interleaved = InherentDiscreteGrid{3}(R, origin; step=step, base=base, unfoldingscheme=:interleaved)

    # Test that coordinate conversions are consistent between schemes
    test_grididx = [(1, 1, 1), (2, 3, 4), (8, 16, 16)]
    for grididx in test_grididx
        coord_fused = QuanticsGrids.grididx_to_origcoord(grid_fused, grididx)
        coord_interleaved = QuanticsGrids.grididx_to_origcoord(grid_interleaved, grididx)
        @test coord_fused == coord_interleaved

        @test QuanticsGrids.origcoord_to_grididx(grid_fused, coord_fused) == grididx
        @test QuanticsGrids.origcoord_to_grididx(grid_interleaved, coord_interleaved) == grididx
    end

    # Test that quantics differ but origcoord results are the same
    grididx = (3, 5, 7)
    quantics_fused = QuanticsGrids.grididx_to_quantics(grid_fused, grididx)
    quantics_interleaved = QuanticsGrids.grididx_to_quantics(grid_interleaved, grididx)

    @test length(quantics_fused) == R
    @test length(quantics_interleaved) == 3 * R
    @test quantics_fused != quantics_interleaved  # Should be different

    # But they should convert to same coordinates
    coord_from_fused = QuanticsGrids.quantics_to_origcoord(grid_fused, quantics_fused)
    coord_from_interleaved = QuanticsGrids.quantics_to_origcoord(grid_interleaved, quantics_interleaved)
    @test coord_from_fused == coord_from_interleaved
end

@testitem "InherentDiscreteGrid high-dimensional grids" begin
    R = 2  # Keep R small for high dimensions to avoid memory issues
    base = 2
    D = 5  # 5-dimensional grid

    origin = ntuple(i -> i, D)  # (1, 2, 3, 4, 5)
    step = ntuple(i -> i, D)    # (1, 2, 3, 4, 5)

    grid = InherentDiscreteGrid{D}(R, origin; step=step, base=base)

    @test ndims(grid) == D
    @test QuanticsGrids.grid_min(grid) == origin
    @test QuanticsGrids.grid_step(grid) == step

    # Test coordinate conversion
    grididx = ntuple(i -> 2, D)  # All 2s
    expected_origcoord = origin .+ step .* (grididx .- 1)
    @test QuanticsGrids.grididx_to_origcoord(grid, grididx) == expected_origcoord
    @test QuanticsGrids.origcoord_to_grididx(grid, expected_origcoord) == grididx

    # Test quantics
    quantics = QuanticsGrids.grididx_to_quantics(grid, grididx)
    @test QuanticsGrids.quantics_to_grididx(grid, quantics) == grididx
    @test length(quantics) == R  # Fused scheme
end

@testitem "InherentDiscreteGrid performance and stress tests" begin
    R = 5
    base = 2
    grid = InherentDiscreteGrid{3}(R, (1, 1, 1); base=base)

    # Test random conversions for consistency
    for _ in 1:100
        # Generate random valid grid indices
        grididx = ntuple(3) do d
            rand(1:base^R)
        end

        # Test roundtrip conversions
        origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
        @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx

        quantics = QuanticsGrids.grididx_to_quantics(grid, grididx)
        @test QuanticsGrids.quantics_to_grididx(grid, quantics) == grididx

        # Test cross-conversions
        origcoord_from_quantics = QuanticsGrids.quantics_to_origcoord(grid, quantics)
        @test origcoord_from_quantics == origcoord

        quantics_from_origcoord = QuanticsGrids.origcoord_to_quantics(grid, origcoord)
        @test quantics_from_origcoord == quantics
    end
end

@testitem "InherentDiscreteGrid grid with zero step (degenerate case)" begin
    # This might be an edge case that should either work or throw a meaningful error
    R = 3
    origin = (5, 10)

    @test_throws AssertionError InherentDiscreteGrid{2}(R, origin; step=(0, 1))
end

@testitem "InherentDiscreteGrid with custom indextable - basic" begin
    # Test with a simple custom indextable
    R = 3
    variablenames = (:x, :y)
    indextable = [
        [(:x, 1)],                    # site 1: x_1
        [(:x, 2)],                    # site 2: x_2
        [(:y, 1)],                    # site 3: y_1
        [(:x, 3), (:y, 2)]           # site 4: x_3, y_2
    ]
    origin = (1, 5)
    step = (2, 3)
    base = 2

    grid = InherentDiscreteGrid(variablenames, indextable;
        origin=origin, step=step, base=base)

    @test grid.variablenames == variablenames
    @test grid.Rs == (3, 2)  # x has 3 quantics, y has 2 quantics
    @test grid.origin == origin
    @test grid.step == step
    @test grid.base == base
    @test length(grid) == 4  # 4 sites in tensor train

    # Test local dimensions: [2^1, 2^1, 2^1, 2^2] = [2, 2, 2, 4]
    @test QuanticsGrids.localdimensions(grid) == [2, 2, 2, 4]

    # Test coordinate conversions work
    grididx = (3, 2)
    origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
    @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx

    # Test quantics conversion
    quantics = QuanticsGrids.grididx_to_quantics(grid, grididx)
    @test length(quantics) == 4
    @test QuanticsGrids.quantics_to_grididx(grid, quantics) == grididx
    @test QuanticsGrids.quantics_to_origcoord(grid, quantics) == origcoord
end

@testitem "InherentDiscreteGrid with complex indextable - mixed sites" begin
    # Test with a more complex indextable with mixed site sizes
    variablenames = (:a, :b, :c)
    indextable = [
        [(:a, 1), (:b, 1)],          # site 1: a_1, b_1 (2 indices)
        [(:c, 1)],                    # site 2: c_1 (1 index)
        [(:a, 2), (:b, 2), (:c, 2)], # site 3: a_2, b_2, c_2 (3 indices)
        [(:a, 3)],                    # site 4: a_3 (1 index)
        [(:b, 3), (:c, 3)]           # site 5: b_3, c_3 (2 indices)
    ]
    origin = (10, 20, 30)
    step = (1, 2, 5)
    base = 3

    grid = InherentDiscreteGrid(variablenames, indextable;
        origin=origin, step=step, base=base)

    @test grid.variablenames == variablenames
    @test grid.Rs == (3, 3, 3)  # Each variable has 3 quantics
    @test grid.origin == origin
    @test grid.step == step
    @test grid.base == base
    @test length(grid) == 5  # 5 sites in tensor train

    # Test local dimensions: [3^2, 3^1, 3^3, 3^1, 3^2] = [9, 3, 27, 3, 9]
    expected_localdims = [9, 3, 27, 3, 9]
    @test QuanticsGrids.localdimensions(grid) == expected_localdims

    # Test that all quantics sites have correct dimensions
    for (i, expected_dim) in enumerate(expected_localdims)
        @test QuanticsGrids.sitedim(grid, i) == expected_dim
    end

    # Test coordinate conversions for various grid indices
    test_grididx = [(1, 1, 1), (2, 3, 1), (9, 27, 27), (5, 10, 15)]
    for grididx in test_grididx
        if all(grididx .<= base .^ grid.Rs)  # Only test valid indices
            origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
            @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx

            quantics = QuanticsGrids.grididx_to_quantics(grid, grididx)
            @test length(quantics) == 5
            @test all(1 .<= quantics[i] .<= expected_localdims[i] for i in 1:5)
            @test QuanticsGrids.quantics_to_grididx(grid, quantics) == grididx
        end
    end
end

@testitem "InherentDiscreteGrid with asymmetric indextable" begin
    # Test with highly asymmetric resolution per variable
    variablenames = (:x, :y, :z)
    indextable = [
        [(:x, 1), (:x, 2), (:x, 3), (:x, 4)],  # site 1: all x quantics (4 indices)
        [(:y, 1)],                               # site 2: y_1 (1 index)
        [(:y, 2)],                               # site 3: y_2 (1 index)
        [(:z, 1), (:z, 2)]                      # site 4: z_1, z_2 (2 indices)
    ]
    origin = (0, 100, -50)
    step = (10, 5, 2)
    base = 2

    grid = InherentDiscreteGrid(variablenames, indextable;
        origin=origin, step=step, base=base)

    @test grid.variablenames == variablenames
    @test grid.Rs == (4, 2, 2)  # x has 4 quantics, y and z have 2 each
    @test length(grid) == 4  # 4 sites

    # Test local dimensions: [2^4, 2^1, 2^1, 2^2] = [16, 2, 2, 4]
    @test QuanticsGrids.localdimensions(grid) == [16, 2, 2, 4]

    # Test that different variables have different resolutions
    max_grididx_x = base^4  # 16
    max_grididx_y = base^2  # 4
    max_grididx_z = base^2  # 4

    # Test boundary cases
    grididx = (max_grididx_x, max_grididx_y, max_grididx_z)
    origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
    expected_origcoord = origin .+ step .* (grididx .- 1)
    @test origcoord == expected_origcoord
    @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx
end

@testitem "InherentDiscreteGrid with single-site indextable" begin
    # Test edge case: all quantics in a single site
    variablenames = (:x, :y)
    indextable = [
        [(:x, 1), (:x, 2), (:x, 3), (:y, 1), (:y, 2)]  # All quantics in one site
    ]
    origin = (5, 10)
    step = (1, 3)
    base = 2

    grid = InherentDiscreteGrid(variablenames, indextable;
        origin=origin, step=step, base=base)

    @test grid.variablenames == variablenames
    @test grid.Rs == (3, 2)  # x has 3 quantics, y has 2 quantics
    @test length(grid) == 1  # Only 1 site

    # Test local dimensions: [2^5] = [32] (5 total quantics in one site)
    @test QuanticsGrids.localdimensions(grid) == [32]

    # Test that conversions still work correctly
    grididx = (4, 3)
    origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
    @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx

    quantics = QuanticsGrids.grididx_to_quantics(grid, grididx)
    @test length(quantics) == 1
    @test 1 <= quantics[1] <= 32
    @test QuanticsGrids.quantics_to_grididx(grid, quantics) == grididx
end

@testitem "InherentDiscreteGrid with maximum fragmentation" begin
    # Test extreme case: each quantics in its own site
    variablenames = (:a, :b)
    indextable = [
        [(:a, 1)],    # site 1: a_1
        [(:a, 2)],    # site 2: a_2
        [(:a, 3)],    # site 3: a_3
        [(:b, 1)],    # site 4: b_1
        [(:b, 2)]     # site 5: b_2
    ]
    origin = (1, 1)
    step = (1, 1)
    base = 3

    grid = InherentDiscreteGrid(variablenames, indextable;
        origin=origin, step=step, base=base)

    @test grid.variablenames == variablenames
    @test grid.Rs == (3, 2)  # a has 3 quantics, b has 2 quantics
    @test length(grid) == 5  # 5 sites

    # Test local dimensions: all sites have single quantics, so [3, 3, 3, 3, 3]
    @test QuanticsGrids.localdimensions(grid) == fill(3, 5)

    # Test conversions
    for grididx in [(1, 1), (3, 3), (9, 9), (27, 9)]
        if all(grididx .<= base .^ grid.Rs)
            origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
            @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx

            quantics = QuanticsGrids.grididx_to_quantics(grid, grididx)
            @test length(quantics) == 5
            @test all(1 .<= quantics .<= 3)
            @test QuanticsGrids.quantics_to_grididx(grid, quantics) == grididx
        end
    end
end

@testitem "InherentDiscreteGrid complex indextable with base != 2" begin
    # Test complex indextable with non-binary base
    variablenames = (:u, :v, :w)
    indextable = [
        [(:u, 1), (:v, 1)],          # site 1: u_1, v_1
        [(:w, 1), (:w, 2)],          # site 2: w_1, w_2
        [(:u, 2)],                    # site 3: u_2
        [(:v, 2), (:u, 3), (:w, 3)]  # site 4: v_2, u_3, w_3
    ]
    origin = (-5, 0, 10)
    step = (3, 7, 2)
    base = 5  # Base 5

    grid = InherentDiscreteGrid(variablenames, indextable;
        origin=origin, step=step, base=base)

    @test grid.variablenames == variablenames
    @test grid.Rs == (3, 2, 3)  # u has 3, v has 2, w has 3 quantics
    @test grid.base == 5
    @test length(grid) == 4

    # Test local dimensions: [5^2, 5^2, 5^1, 5^3] = [25, 25, 5, 125]
    @test QuanticsGrids.localdimensions(grid) == [25, 25, 5, 125]

    # Test with base-5 specific values
    grididx = (1, 1, 1)  # Minimum
    origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
    @test origcoord == origin
    @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx

    # Test maximum valid indices
    max_grididx = (base^3, base^2, base^3)  # (125, 25, 125)
    max_origcoord = QuanticsGrids.grididx_to_origcoord(grid, max_grididx)
    expected_max = origin .+ step .* (max_grididx .- 1)
    @test max_origcoord == expected_max
    @test QuanticsGrids.origcoord_to_grididx(grid, max_origcoord) == max_grididx

    # Test quantics conversions
    test_grididx = (10, 5, 20)
    quantics = QuanticsGrids.grididx_to_quantics(grid, test_grididx)
    @test length(quantics) == 4
    @test all(1 .<= quantics[1] .<= 25)   # site 1
    @test all(1 .<= quantics[2] .<= 25)   # site 2
    @test all(1 .<= quantics[3] .<= 5)    # site 3
    @test all(1 .<= quantics[4] .<= 125)  # site 4
    @test QuanticsGrids.quantics_to_grididx(grid, quantics) == test_grididx
end

@testitem "InherentDiscreteGrid comprehensive indextable test" begin
    # Test a realistic complex case similar to what might be used in practice
    variablenames = (:k_x, :k_y, :ω, :t)
    indextable = [
        [(:k_x, 1), (:k_y, 1)],           # site 1: spatial momenta first bits
        [(:ω, 1), (:t, 1)],               # site 2: frequency and time first bits
        [(:k_x, 2)],                      # site 3: k_x second bit
        [(:k_y, 2), (:ω, 2)],             # site 4: k_y and ω second bits
        [(:t, 2), (:t, 3)],               # site 5: time higher resolution
        [(:k_x, 3), (:k_y, 3), (:ω, 3)]  # site 6: spatial and frequency high bits
    ]
    origin = (0, 0, -10, 0)
    step = (1, 1, 2, 5)
    base = 2

    grid = InherentDiscreteGrid(variablenames, indextable;
        origin=origin, step=step, base=base)

    @test grid.variablenames == variablenames
    @test grid.Rs == (3, 3, 3, 3)  # All variables have 3 quantics
    @test length(grid) == 6  # 6 sites

    # Test local dimensions: [2^2, 2^2, 2^1, 2^2, 2^2, 2^3] = [4, 4, 2, 4, 4, 8]
    expected_localdims = [4, 4, 2, 4, 4, 8]
    @test QuanticsGrids.localdimensions(grid) == expected_localdims

    # Test random conversions for consistency
    for _ in 1:50
        # Generate random valid grid indices
        grididx = ntuple(4) do d
            rand(1:base^3)  # Each dimension has 3 quantics, so max index is 2^3 = 8
        end

        # Test all conversion roundtrips
        origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
        @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx

        quantics = QuanticsGrids.grididx_to_quantics(grid, grididx)
        @test length(quantics) == 6
        @test all(1 .<= quantics[i] .<= expected_localdims[i] for i in 1:6)
        @test QuanticsGrids.quantics_to_grididx(grid, quantics) == grididx

        # Cross conversions
        @test QuanticsGrids.quantics_to_origcoord(grid, quantics) == origcoord
        @test QuanticsGrids.origcoord_to_quantics(grid, origcoord) == quantics
    end
end

@testitem "InherentDiscreteGrid indextable validation" begin
    # Test that constructor properly validates indextable consistency
    variablenames = (:x, :y)

    # Valid indextable
    valid_indextable = [
        [(:x, 1), (:y, 1)],
        [(:x, 2)],
        [(:y, 2)]
    ]

    grid = InherentDiscreteGrid(variablenames, valid_indextable;
        origin=(1, 1), step=(1, 1))
    @test grid.Rs == (2, 2)  # x and y each have 2 quantics

    # Test that unknown variables are detected
    invalid_indextable_unknown = [
        [(:x, 1), (:y, 1)],
        [(:z, 1)]  # z is not in variablenames
    ]

    @test_throws AssertionError InherentDiscreteGrid(
        variablenames, invalid_indextable_unknown;
        origin=(1, 1), step=(1, 1))

    invalid_indextable_missing_index = [
        [(:x, 1), (:y, 1)],
        [(:x, 3), (:y, 2)]
    ]

    @test_throws AssertionError InherentDiscreteGrid(
        variablenames, invalid_indextable_missing_index)

    invalid_indextable_repeated_index = [
        [(:x, 1), (:y, 1)],
        [(:x, 1), (:y, 2)]
    ]

    @test_throws AssertionError InherentDiscreteGrid(
        variablenames, invalid_indextable_repeated_index)
end

@testitem "InherentDiscreteGrid constructor with variablenames and Rs" begin
    # Test basic constructor with explicit variable names and resolutions
    variablenames = (:k_x, :k_y, :ω)
    Rs = (4, 3, 5)
    origin = (1, 10, -5)
    step = (2, 1, 3)

    # Test with default parameters
    grid1 = InherentDiscreteGrid(variablenames, Rs; origin=origin, step=step)
    @test grid1.variablenames == variablenames
    @test grid1.Rs == Rs
    @test grid1.base == 2  # default base
    @test grid1.origin == origin
    @test grid1.step == step
    @test ndims(grid1) == 3

    # Test with custom parameters
    base = 3
    unfoldingscheme = :interleaved

    grid2 = InherentDiscreteGrid(variablenames, Rs;
        origin=origin,
        step=step,
        base=base,
        unfoldingscheme=unfoldingscheme)

    @test grid2.variablenames == variablenames
    @test grid2.Rs == Rs
    @test grid2.base == base
    @test grid2.origin == origin
    @test grid2.step == step

    # Test that functionality works
    test_grididx = (2, 3, 4)
    origcoord = QuanticsGrids.grididx_to_origcoord(grid2, test_grididx)
    back_to_grididx = QuanticsGrids.origcoord_to_grididx(grid2, origcoord)
    @test back_to_grididx == test_grididx

    # Test quantics conversion
    quantics = QuanticsGrids.grididx_to_quantics(grid2, test_grididx)
    @test QuanticsGrids.quantics_to_grididx(grid2, quantics) == test_grididx
end

@testitem "InherentDiscreteGrid constructor with Rs tuple only" begin
    # Test constructor that auto-generates variable names
    Rs = (3, 4, 2)
    origin = (5, 0, -10)
    step = (1, 2, 5)
    base = 3

    grid = InherentDiscreteGrid(Rs; origin=origin, step=step, base=base)

    @test grid.Rs == Rs
    @test grid.origin == origin
    @test grid.step == step
    @test grid.base == base
    @test ndims(grid) == 3

    # Should auto-generate variable names as symbols
    expected_names = (Symbol(1), Symbol(2), Symbol(3))
    @test grid.variablenames == expected_names

    # Test that grid works correctly
    grididx = (2, 3, 1)
    origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
    @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx
end

@testitem "InherentDiscreteGrid 1D constructor with R::Int" begin
    # Test 1D constructor similar to DiscretizedGrid(R, lower, upper)
    R = 5
    origin = 10
    step = 3
    base = 2

    grid = InherentDiscreteGrid(R, origin; step=step, base=base)

    @test grid.Rs == (R,)
    @test grid.origin == (origin,)
    @test grid.step == (step,)
    @test grid.base == base
    @test ndims(grid) == 1

    # Should auto-generate 1D variable name
    @test grid.variablenames == (Symbol(1),)

    # Test conversions work for 1D case (should return scalars)
    grididx = 3
    origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
    @test origcoord isa Int
    @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx

    # Test expected coordinate calculation
    expected_coord = origin + step * (grididx - 1)
    @test origcoord == expected_coord
end

@testitem "InherentDiscreteGrid parametric constructor InherentDiscreteGrid{D}" begin
    # Test parametric type constructor like DiscretizedGrid{D}
    D = 4
    R = 3
    origin = ntuple(i -> i * 5, D)  # (5, 10, 15, 20)
    step = ntuple(i -> i, D)        # (1, 2, 3, 4)
    base = 3

    grid = InherentDiscreteGrid{D}(R, origin; step=step, base=base, unfoldingscheme=:interleaved)

    @test ndims(grid) == D
    @test grid.Rs == ntuple(_ -> R, D)  # All dimensions have same R
    @test grid.origin == origin
    @test grid.step == step
    @test grid.base == base

    # Should auto-generate variable names
    expected_names = ntuple(Symbol, D)
    @test grid.variablenames == expected_names

    # Test that grid works for high dimensions
    grididx = ntuple(i -> 2, D)
    origcoord = QuanticsGrids.grididx_to_origcoord(grid, grididx)
    @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == grididx
end

@testitem "InherentDiscreteGrid default parameter behavior" begin
    # Test various combinations of default parameters
    Rs = (3, 2)

    # Test with minimal parameters (all defaults)
    grid_minimal = InherentDiscreteGrid(Rs)
    @test grid_minimal.Rs == Rs
    @test grid_minimal.base == 2  # default
    @test grid_minimal.origin == (1, 1)  # default origin should be (1, 1, ...)
    @test grid_minimal.step == (1, 1)    # default step should be (1, 1, ...)

    # Test with partial parameters
    custom_origin = (5, -2)
    grid_partial = InherentDiscreteGrid(Rs; origin=custom_origin)
    @test grid_partial.origin == custom_origin
    @test grid_partial.step == (1, 1)  # should still be default
    @test grid_partial.base == 2       # should still be default

    # Test step can be specified as single value or tuple
    grid_step_single = InherentDiscreteGrid(Rs; step=3)
    @test grid_step_single.step == (3, 3)  # should broadcast to all dimensions

    grid_step_tuple = InherentDiscreteGrid(Rs; step=(2, 5))
    @test grid_step_tuple.step == (2, 5)

    # Test origin can be specified as single value or tuple
    grid_origin_single = InherentDiscreteGrid(Rs; origin=7)
    @test grid_origin_single.origin == (7, 7)  # should broadcast

    grid_origin_tuple = InherentDiscreteGrid(Rs; origin=(-1, 3))
    @test grid_origin_tuple.origin == (-1, 3)
end

@testitem "InherentDiscreteGrid constructor error handling" begin
    # Test that constructors properly validate input parameters

    # Test invalid base
    @test_throws AssertionError InherentDiscreteGrid((3, 2); base=1)  # base must be > 1
    @test_throws AssertionError InherentDiscreteGrid((3, 2); base=0)

    # Test mismatched dimensions
    @test_throws MethodError InherentDiscreteGrid((3, 2); origin=(1, 1, 1))  # wrong origin length
    @test_throws MethodError InherentDiscreteGrid((3, 2); step=(1, 1, 1))    # wrong step length

    # Test negative R values
    @test_throws AssertionError InherentDiscreteGrid((-1, 2))  # negative R
    @test_throws AssertionError InherentDiscreteGrid((3, -2))

    # Test invalid unfoldingscheme
    @test_throws AssertionError InherentDiscreteGrid((3, 2); unfoldingscheme=:invalid)

    # Test zero R (might be edge case depending on implementation)
    try
        grid_zero = InherentDiscreteGrid((0, 2))
        # If it works, verify it behaves sensibly
        @test grid_zero.Rs == (0, 2)
    catch err
        # If it throws an error, verify it's meaningful
        @test err isa AssertionError || err isa DomainError
    end
end

@testitem "InherentDiscreteGrid constructor compatibility with existing patterns" begin
    # Test that constructors work with patterns used in existing codebase

    # Pattern 1: Like InherentDiscreteGrid{d}(R, origin; kwargs...)
    R = 4
    origin_2d = (1, 5)
    step_2d = (2, 3)

    grid_pattern1 = InherentDiscreteGrid{2}(R, origin_2d; step=step_2d, base=3)
    @test ndims(grid_pattern1) == 2
    @test grid_pattern1.Rs == (R, R)
    @test grid_pattern1.origin == origin_2d
    @test grid_pattern1.step == step_2d

    # Pattern 2: Like DiscretizedGrid((R1, R2, ...); kwargs...)
    Rs_pattern2 = (3, 5, 2)
    grid_pattern2 = InherentDiscreteGrid(Rs_pattern2; origin=(0, 0, 0), unfoldingscheme=:interleaved)
    @test grid_pattern2.Rs == Rs_pattern2

    # Pattern 3: Single R for 1D
    R_1d = 6
    origin_1d = 10
    grid_pattern3 = InherentDiscreteGrid(R_1d, origin_1d; step=2, base=4)
    @test ndims(grid_pattern3) == 1
    @test grid_pattern3.Rs == (R_1d,)
    @test grid_pattern3.origin == (origin_1d,)
    @test grid_pattern3.base == 4

    # Test that all patterns produce working grids
    for grid in [grid_pattern1, grid_pattern2, grid_pattern3]
        # Test basic coordinate conversion works
        if ndims(grid) == 1
            test_idx = 2
        elseif ndims(grid) == 2
            test_idx = (2, 3)
        else
            test_idx = (2, 3, 1)
        end

        origcoord = QuanticsGrids.grididx_to_origcoord(grid, test_idx)
        @test QuanticsGrids.origcoord_to_grididx(grid, origcoord) == test_idx
    end
end

@testitem "InherentDiscreteGrid constructor with mixed parameter types" begin
    # Test constructors handle mixed integer/tuple parameters correctly

    # Single R with multi-dimensional origin/step
    R = 3
    origin_multi = (5, 10, 15)
    step_multi = (1, 2, 3)

    grid_mixed = InherentDiscreteGrid{3}(R, origin_multi; step=step_multi)
    @test grid_mixed.Rs == (R, R, R)  # R should be replicated
    @test grid_mixed.origin == origin_multi
    @test grid_mixed.step == step_multi

    # Multi-R with single origin/step values
    Rs_multi = (2, 4, 3)
    origin_single = 7
    step_single = 2

    grid_broadcast = InherentDiscreteGrid(Rs_multi; origin=origin_single, step=step_single)
    @test grid_broadcast.Rs == Rs_multi
    @test grid_broadcast.origin == (origin_single, origin_single, origin_single)
    @test grid_broadcast.step == (step_single, step_single, step_single)

    # Test functionality works correctly
    grididx = (1, 2, 3)
    expected_origcoord = (origin_single, origin_single + step_single * 1, origin_single + step_single * 2)
    actual_origcoord = QuanticsGrids.grididx_to_origcoord(grid_broadcast, grididx)
    @test actual_origcoord == expected_origcoord
end
