@testitem "localdimensions" begin
    R = 4
    # Basic test with default base=2
    grid = NewDiscretizedGrid{1}(R; lower_bound=(0.0,), upper_bound=(1.0,))
    @test QuanticsGrids.localdimensions(grid) == fill(2, R)

    grid_2d = NewDiscretizedGrid((R, R); lower_bound=(0.0, 0.0), upper_bound=(1.0, 1.0))
    @test QuanticsGrids.localdimensions(grid_2d) == fill(2, 2R)

    # Test with different base values
    grid_base3 = NewDiscretizedGrid{1}(R; lower_bound=(0.0,), upper_bound=(1.0,), base=3)
    @test QuanticsGrids.localdimensions(grid_base3) == fill(3, R)

    grid_base4 = NewDiscretizedGrid{1}(R; lower_bound=(0.0,), upper_bound=(1.0,), base=4)
    @test QuanticsGrids.localdimensions(grid_base4) == fill(4, R)

    # Test with custom indextable - sites with multiple quantics indices
    variablenames = (:x, :y)
    # Create indextable where some sites have multiple quantics indices
    indextable = [
        [(:x, 1), (:y, 1)],  # site 1: 2 indices
        [(:x, 2)],           # site 2: 1 index  
        [(:y, 2)],           # site 3: 1 index
        [(:x, 3), (:y, 3), (:x, 4)]  # site 4: 3 indices
    ]

    # Test with base=2: Quan  ticsGrids.localdimensions should be [2^2, 2^1, 2^1, 2^3] = [4, 2, 2, 8]
    grid_custom_base2 = NewDiscretizedGrid(variablenames, indextable;
        lower_bound=(0.0, 0.0), upper_bound=(1.0, 1.0), base=2)
    @test QuanticsGrids.localdimensions(grid_custom_base2) == [4, 2, 2, 8]

    # Test with base=3: Quan  ticsGrids.localdimensions should be [3^2, 3^1, 3^1, 3^3] = [9, 3, 3, 27]
    grid_custom_base3 = NewDiscretizedGrid(variablenames, indextable;
        lower_bound=(0.0, 0.0), upper_bound=(1.0, 1.0), base=3)
    @test QuanticsGrids.localdimensions(grid_custom_base3) == [9, 3, 3, 27]

    # Test edge case: single site with multiple indices
    single_site_indextable = [
        [(:x, 1), (:x, 2), (:x, 3)]
    ]
    grid_single_site = NewDiscretizedGrid((:x,), single_site_indextable;
        lower_bound=(0.0,), upper_bound=(1.0,), base=2)
    @test QuanticsGrids.localdimensions(grid_single_site) == [8]  # 2^3 = 8

    # Test with mixed site sizes in higher dimension
    complex_indextable = [
        [(:x, 1)],                          # site 1: 1 index -> base^1
        [(:y, 1), (:z, 1)],                 # site 2: 2 indices -> base^2
        [(:x, 2), (:y, 2), (:z, 2)],        # site 3: 3 indices -> base^3
        [(:x, 3), (:y, 3)]                  # site 4: 2 indices -> base^2
    ]
    grid_complex = NewDiscretizedGrid((:x, :y, :z), complex_indextable;
        lower_bound=(0.0, 0.0, 0.0), upper_bound=(1.0, 1.0, 1.0), base=2)
    @test QuanticsGrids.localdimensions(grid_complex) == [2, 4, 8, 4]  # [2^1, 2^2, 2^3, 2^2]
end

@testitem "quanticsfunction" begin
    R = 4
    grid = NewDiscretizedGrid{1}(R; lower_bound=(0.0,), upper_bound=(1.0,))

    # Test that quanticsfunction should work with NewDiscretizedGrid
    fx(x) = exp(-x)
    fq = QuanticsGrids.quanticsfunction(Float64, grid, fx)

    # Test basic functionality - should convert quantics to coordinates and apply function
    @test fq(fill(1, R)) == fx(0)  # quantics [1,1,1,1] should map to x=0.0
    @test fq(fill(2, R)) == fx(1 - 1 / 2^R)  # quantics [2,2,2,2] should map to near x=1.0

    # Test 2D case
    grid_2d = NewDiscretizedGrid((R, R); lower_bound=(0.0, 0.0), upper_bound=(1.0, 1.0))
    f2d(x, y) = x + y
    fq_2d = QuanticsGrids.quanticsfunction(Float64, grid_2d, f2d)

    # Test 2D functionality
    quantics_origin = fill(1, length(grid_2d))  # Should map to (0.0, 0.0)
    @test fq_2d(quantics_origin) == f2d(0.0, 0.0)

    # Test with custom indextable
    variablenames = (:x, :y)
    indextable = [
        [(:x, 1)],
        [(:x, 2)],
        [(:y, 1)],
        [(:y, 2)]
    ]
    grid_custom = NewDiscretizedGrid(variablenames, indextable;
        lower_bound=(0.0, 0.0), upper_bound=(1.0, 1.0))
    fq_custom = QuanticsGrids.quanticsfunction(Float64, grid_custom, f2d)
    @test fq_custom([1, 1, 1, 1]) == f2d(0.0, 0.0)

    # CHALLENGING TEST: Complex multi-variable function with mixed site sizes and base=3
    # Create a complex indextable with different numbers of quantics per site
    complex_variablenames = (:x, :y, :z)
    complex_indextable = [
        [(:x, 1)],
        [(:y, 1), (:z, 1)],
        [(:x, 2), (:y, 2)],
        [(:z, 2), (:x, 3), (:y, 3)]
    ]

    # Create grid with base=3, non-unit domain
    grid_complex = NewDiscretizedGrid(complex_variablenames, complex_indextable;
        lower_bound=(-1.0, 2.0, 0.5), upper_bound=(3.0, 5.0, 2.0), base=3)

    # Define a challenging multi-variable function with cross-terms and oscillations
    f_complex(x, y, z) = sin(π * x) * cos(2π * y) * exp(-z) + x^2 * y - z^3 + x * y * z
    fq_complex = QuanticsGrids.quanticsfunction(Float64, grid_complex, f_complex)

    # Test corner cases: minimum coordinates (all quantics = 1)
    @test fq_complex([1, 1, 1, 1]) ≈ f_complex(-1.0, 2.0, 0.5)

    # Test maximum coordinates (all quantics at max value for each site's local dimension)
    @test fq_complex([3, 9, 9, 27]) ≈ f_complex(QuanticsGrids.grid_max(grid_complex)...)

    # Test intermediate values with specific quantics combinations
    # Choose quantics that should map to approximately middle values
    quantics_mid = [2, 5, 5, 14]  # Should map to roughly middle of each dimension
    coord_mid = quantics_to_origcoord(grid_complex, quantics_mid)
    @test fq_complex(quantics_mid) ≈ f_complex(coord_mid...)
end

@testitem "unfoldingscheme" begin
    R = 4

    # NewDiscretizedGrid should support unfoldingscheme parameter and property
    grid_fused = NewDiscretizedGrid{1}(R; lower_bound=(0.0,), upper_bound=(1.0,), unfoldingscheme=:fused)
    grid_interleaved = NewDiscretizedGrid{1}(R; lower_bound=(0.0,), upper_bound=(1.0,), unfoldingscheme=:interleaved)

    # Test 2D grids with different unfolding schemes
    grid_2d_fused = NewDiscretizedGrid((R, R); lower_bound=(0.0, 0.0), upper_bound=(1.0, 1.0), unfoldingscheme=:fused)
    grid_2d_interleaved = NewDiscretizedGrid((R, R); lower_bound=(0.0, 0.0), upper_bound=(1.0, 1.0), unfoldingscheme=:interleaved)

    # Test that unfoldingscheme affects quantics behavior
    # For 2D grids: fused should have length R, interleaved should have length 2*R
    grididx = (2, 2)
    quantics_fused = grididx_to_quantics(grid_2d_fused, grididx)
    quantics_interleaved = grididx_to_quantics(grid_2d_interleaved, grididx)

    @test length(quantics_fused) == R  # fused: length R
    @test length(quantics_interleaved) == 2 * R  # interleaved: length 2*R for 2D

    # Test that both produce same coordinate when converted back
    coord_fused = quantics_to_origcoord(grid_2d_fused, quantics_fused)
    coord_interleaved = quantics_to_origcoord(grid_2d_interleaved, quantics_interleaved)
    @test coord_fused == coord_interleaved

    # Test default unfoldingscheme (should be :fused)
    grid_default = NewDiscretizedGrid{1}(R; lower_bound=(0.0,), upper_bound=(1.0,))

    # Test localdimensions with different unfolding schemes
    @test QuanticsGrids.localdimensions(grid_2d_fused) == fill(4, R)  # 2^2 = 4 for fused
    @test QuanticsGrids.localdimensions(grid_2d_interleaved) == fill(2, 2 * R)  # 2 for interleaved
end

@testitem "includeendpoint" begin
    R = 4

    # NewDiscretizedGrid should support includeendpoint parameter and property
    grid_without = NewDiscretizedGrid{1}(R; lower_bound=(0.0,), upper_bound=(1.0,), includeendpoint=false)
    grid_with = NewDiscretizedGrid{1}(R; lower_bound=(0.0,), upper_bound=(1.0,), includeendpoint=true)

    # Test that includeendpoint affects grid behavior
    @test QuanticsGrids.grid_max(grid_without) ≈ 1.0 - (1.0 / 2^R)
    @test QuanticsGrids.grid_max(grid_with) == 1.0

    # Test grid_step calculation differences
    step_without = QuanticsGrids.grid_step(grid_without)
    step_with = QuanticsGrids.grid_step(grid_with)

    @test step_without ≈ 1.0 / 2^R
    @test step_with ≈ 1.0 / (2^R - 1)

    # Test that includeendpoint=true allows accessing the upper bound
    @test grididx_to_origcoord(grid_with, 2^R) == 1.0
    @test grididx_to_origcoord(grid_without, 2^R) ≈ 1.0 - (1.0 / 2^R)

    # Test 2D case
    grid_2d_without = NewDiscretizedGrid((R, R);
        lower_bound=(0.0, 0.0),
        upper_bound=(2.0, 3.0),
        includeendpoint=false)
    grid_2d_with = NewDiscretizedGrid((R, R);
        lower_bound=(0.0, 0.0),
        upper_bound=(2.0, 3.0),
        includeendpoint=true)

    # Test grid_max behavior in 2D
    max_without = QuanticsGrids.grid_max(grid_2d_without)
    max_with = QuanticsGrids.grid_max(grid_2d_with)

    @test max_without == (2.0 - 2.0 / 2^R, 3.0 - 3.0 / 2^R)
    @test max_with == (2.0, 3.0)

    # Test that the upper bound point is accessible when includeendpoint=true
    @test grididx_to_origcoord(grid_2d_with, (2^R, 2^R)) == (2.0, 3.0)

    # CHALLENGING TEST: Complex indextable with mixed site sizes and includeendpoint behavior
    # This tests how includeendpoint interacts with custom indextables having different local dimensions
    complex_variablenames = (:x, :y, :z)
    complex_indextable = [
        [(:x, 1)],
        [(:y, 1), (:z, 1)],
        [(:x, 2), (:y, 2), (:z, 2)],
        [(:x, 3), (:y, 3)]
    ]

    # Test with base=3 and asymmetric domain bounds
    lower_bound = (-2.5, 1.2, 0.0)
    upper_bound = (4.7, 8.1, 3.3)

    grid_complex_without = NewDiscretizedGrid(complex_variablenames, complex_indextable;
        lower_bound, upper_bound, base=3, includeendpoint=false)
    grid_complex_with = NewDiscretizedGrid(complex_variablenames, complex_indextable;
        lower_bound, upper_bound, base=3, includeendpoint=true)

    # Test that local dimensions are same regardless of includeendpoint
    local_dims_without = QuanticsGrids.localdimensions(grid_complex_without)
    local_dims_with = QuanticsGrids.localdimensions(grid_complex_with)
    expected_local_dims = [3, 9, 27, 9]  # [3^1, 3^2, 3^3, 3^2]
    @test local_dims_without == expected_local_dims
    @test local_dims_with == expected_local_dims

    # Test grid_max behavior with complex grids
    max_without = QuanticsGrids.grid_max(grid_complex_without)
    max_with = QuanticsGrids.grid_max(grid_complex_with)

    # Without includeendpoint: should be slightly less than upper_bound bound
    expected_max_without = (
        upper_bound[1] - (upper_bound[1] - lower_bound[1]) / 3^3,  # x dimension uses 3 quantics
        upper_bound[2] - (upper_bound[2] - lower_bound[2]) / 3^3,  # y dimension uses 3 quantics  
        upper_bound[3] - (upper_bound[3] - lower_bound[3]) / 3^2   # z dimension uses 2 quantics
    )
    @test all(max_without .≈ expected_max_without)
    @test max_with == upper_bound

    # Test maximum quantics combinations for each site's local dimension
    max_quantics_without = [3, 9, 27, 9]  # maximum for each site
    max_quantics_with = [3, 9, 27, 9]     # same maximum values

    coord_max_without = quantics_to_origcoord(grid_complex_without, max_quantics_without)
    coord_max_with = quantics_to_origcoord(grid_complex_with, max_quantics_with)

    # With includeendpoint=true, max quantics should reach exact upper_bound bound
    @test all(coord_max_with .≈ upper_bound)

    # With includeendpoint=false, max quantics should be slightly less than upper_bound bound
    @test all(coord_max_without .< upper_bound)
    @test all(coord_max_without .≈ expected_max_without)

    # Test quanticsfunction behavior with boundary values
    # Define a function that's sensitive to boundary effects
    boundary_sensitive_func(x, y, z) = (x - upper_bound[1])^4 + (y - upper_bound[2])^4 + (z - upper_bound[3])^4

    fq_without = QuanticsGrids.quanticsfunction(Float64, grid_complex_without, boundary_sensitive_func)
    fq_with = QuanticsGrids.quanticsfunction(Float64, grid_complex_with, boundary_sensitive_func)

    # Test that includeendpoint=true allows exact evaluation at upper_bound bound
    @test fq_with(max_quantics_with) ≈ boundary_sensitive_func(upper_bound...)

    # Test that includeendpoint=false gives a non-zero value (since it can't reach exact upper_bound bound)
    result_without = fq_without(max_quantics_without)
    @test result_without > 0  # Should be positive since we can't reach the upper_bound bound exactly
    @test result_without ≈ boundary_sensitive_func(coord_max_without...)

    # Test step size consistency across different sites
    # For complex indextables, step sizes should vary by dimension based on quantics count
    step_without = QuanticsGrids.grid_step(grid_complex_without)
    step_with = QuanticsGrids.grid_step(grid_complex_with)

    # Expected step sizes for each dimension
    expected_step_without = (
        (upper_bound[1] - lower_bound[1]) / 3^3,  # x: uses 3 quantics per site
        (upper_bound[2] - lower_bound[2]) / 3^3,  # y: uses 3 quantics per site
        (upper_bound[3] - lower_bound[3]) / 3^2   # z: uses 2 quantics per site
    )
    expected_step_with = (
        (upper_bound[1] - lower_bound[1]) / (3^3 - 1),  # x: adjusted for includeendpoint
        (upper_bound[2] - lower_bound[2]) / (3^3 - 1),  # y: adjusted for includeendpoint
        (upper_bound[3] - lower_bound[3]) / (3^2 - 1)   # z: adjusted for includeendpoint
    )

    @test all(step_without .≈ expected_step_without)
    @test all(step_with .≈ expected_step_with)

    # Test that grid_origin is unaffected by includeendpoint (should always be lower_bound)
    @test QuanticsGrids.grid_origin(grid_complex_without) == lower_bound
    @test QuanticsGrids.grid_origin(grid_complex_with) == lower_bound

    # Test quantics-to-coordinate conversion consistency at boundaries
    # Minimum quantics (all 1s) should always map to lower_bound regardless of includeendpoint
    min_quantics = [1, 1, 1, 1]
    @test all(quantics_to_origcoord(grid_complex_without, min_quantics) .≈ lower_bound)
    @test all(quantics_to_origcoord(grid_complex_with, min_quantics) .≈ lower_bound)

    # Test intermediate quantics values behave consistently
    mid_quantics = [2, 5, 14, 5]  # roughly middle values for each site
    coord_mid_without = quantics_to_origcoord(grid_complex_without, mid_quantics)
    coord_mid_with = quantics_to_origcoord(grid_complex_with, mid_quantics)

    # Coordinates should be different due to different step sizes
    @test coord_mid_without != coord_mid_with

    # But both should be within the domain bounds
    @test all(lower_bound .<= coord_mid_without .<= upper_bound)
    @test all(lower_bound .<= coord_mid_with .<= upper_bound)
end

@testitem "boundary error handling" begin
    R = 5
    g = NewDiscretizedGrid{1}(R)
    @test_throws AssertionError QuanticsGrids.origcoord_to_grididx(g, -0.1)
    @test_throws AssertionError QuanticsGrids.origcoord_to_grididx(g, 1.1)
end

@testitem "large grids" begin
    R = 62
    g = NewDiscretizedGrid{1}(R)
    @test grididx_to_quantics(g, 2^R) == fill(2, R)
end

@testitem "grid representation conversion _NEW_" begin
    R = 10
    reprs = [:grididx, :quantics, :origcoord]

    varnames = (:x, :y)
    tupletable = [
        [(:x, 3), (:y, 2)],
        [(:y, 1), (:x, 1)],
        [(:x, 2)],
    ]

    testset = [
        (NewDiscretizedGrid{1}(R, -3.2, 4.8), 2),
        (NewDiscretizedGrid{2}(R, (-2.3, 1.2), (9.0, 3.5)), (2, 3)),
        (NewDiscretizedGrid(varnames, tupletable), (7, 2)),
    ]

    for (grid, ini_grididx) in testset

        data = Dict{Symbol,Any}()
        transforms = Dict(
            (:grididx, :quantics) => grididx_to_quantics,
            (:grididx, :origcoord) => grididx_to_origcoord,
            (:quantics, :grididx) => quantics_to_grididx,
            (:quantics, :origcoord) => quantics_to_origcoord,
            (:origcoord, :grididx) => QuanticsGrids.origcoord_to_grididx,
            (:origcoord, :quantics) => origcoord_to_quantics,
        )

        data[:grididx] = ini_grididx
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
                    data[dst] = @inferred transforms[(src, dst)](grid, data[src])
                    newdata = true
                end
            end
            if newdata == false
                break
            end
        end
    end
end

@testitem "step size and origin" begin
    R = 4
    g = NewDiscretizedGrid{1}(R)
    @test QuanticsGrids.grid_step(g) == 1.0 / 2^R
    @test QuanticsGrids.grid_origin(g) == 0.0
end

@testitem "grid_origcoords" begin
    # Test 1D grid with default parameters
    R = 4
    grid_1d = NewDiscretizedGrid{1}(R; lower_bound=(0.0,), upper_bound=(1.0,))
    coords_1d = QuanticsGrids.grid_origcoords(grid_1d, 1)

    # Should return a range with correct length and bounds
    @test length(coords_1d) == 2^R
    @test coords_1d[1] ≈ QuanticsGrids.grid_min(grid_1d)
    @test coords_1d[end] ≈ QuanticsGrids.grid_max(grid_1d)

    # Test that coordinates match grididx_to_origcoord for all indices
    for i in 1:2^R
        @test coords_1d[i] ≈ grididx_to_origcoord(grid_1d, i)
    end

    # Test 2D grid
    grid_2d = NewDiscretizedGrid((3, 4); lower_bound=(-1.0, 2.0), upper_bound=(5.0, 8.0))

    # Test first dimension
    coords_x = QuanticsGrids.grid_origcoords(grid_2d, 1)
    @test length(coords_x) == 2^3
    @test coords_x[1] ≈ QuanticsGrids.grid_min(grid_2d)[1]
    @test coords_x[end] ≈ QuanticsGrids.grid_max(grid_2d)[1]

    # Test second dimension  
    coords_y = QuanticsGrids.grid_origcoords(grid_2d, 2)
    @test length(coords_y) == 2^4
    @test coords_y[1] ≈ QuanticsGrids.grid_min(grid_2d)[2]
    @test coords_y[end] ≈ QuanticsGrids.grid_max(grid_2d)[2]

    # Verify coordinates match grididx_to_origcoord
    for i in 1:2^3
        @test coords_x[i] ≈ grididx_to_origcoord(grid_2d, (i, 1))[1]
    end
    for j in 1:2^4
        @test coords_y[j] ≈ grididx_to_origcoord(grid_2d, (1, j))[2]
    end

    # Test with different base
    grid_base3 = NewDiscretizedGrid{1}(3; lower_bound=(0.0,), upper_bound=(2.7,), base=3)
    coords_base3 = QuanticsGrids.grid_origcoords(grid_base3, 1)
    @test length(coords_base3) == 3^3
    @test coords_base3[1] ≈ QuanticsGrids.grid_min(grid_base3)
    @test coords_base3[end] ≈ QuanticsGrids.grid_max(grid_base3)

    # Test with includeendpoint=true
    grid_endpoint = NewDiscretizedGrid{1}(3; lower_bound=(0.0,), upper_bound=(1.0,), includeendpoint=true)
    coords_endpoint = QuanticsGrids.grid_origcoords(grid_endpoint, 1)
    @test length(coords_endpoint) == 2^3
    @test coords_endpoint[1] ≈ 0.0
    @test coords_endpoint[end] ≈ 1.0  # Should reach exact upper bound

    # Test with custom indextable
    variablenames = (:x, :y, :z)
    indextable = [
        [(:x, 1), (:y, 1)],    # site 1: 2 indices
        [(:z, 1)],             # site 2: 1 index
        [(:x, 2), (:y, 2), (:z, 2)],  # site 3: 3 indices
        [(:x, 3), (:y, 3)]     # site 4: 2 indices
    ]
    grid_custom = NewDiscretizedGrid(variablenames, indextable;
        lower_bound=(-2.0, 1.0, 0.5), upper_bound=(3.0, 4.0, 2.0), base=2)

    # Test each dimension
    coords_x_custom = QuanticsGrids.grid_origcoords(grid_custom, 1)  # x dimension (3 quantics)
    coords_y_custom = QuanticsGrids.grid_origcoords(grid_custom, 2)  # y dimension (3 quantics)
    coords_z_custom = QuanticsGrids.grid_origcoords(grid_custom, 3)  # z dimension (2 quantics)

    @test length(coords_x_custom) == 2^3  # x has 3 quantics
    @test length(coords_y_custom) == 2^3  # y has 3 quantics
    @test length(coords_z_custom) == 2^2  # z has 2 quantics

    # Test bounds
    @test coords_x_custom[1] ≈ QuanticsGrids.grid_min(grid_custom)[1]
    @test coords_y_custom[1] ≈ QuanticsGrids.grid_min(grid_custom)[2]
    @test coords_z_custom[1] ≈ QuanticsGrids.grid_min(grid_custom)[3]
    @test coords_x_custom[end] ≈ QuanticsGrids.grid_max(grid_custom)[1]
    @test coords_y_custom[end] ≈ QuanticsGrids.grid_max(grid_custom)[2]
    @test coords_z_custom[end] ≈ QuanticsGrids.grid_max(grid_custom)[3]

    # Test error handling - dimension out of bounds
    @test_throws AssertionError QuanticsGrids.grid_origcoords(grid_1d, 0)
    @test_throws AssertionError QuanticsGrids.grid_origcoords(grid_1d, 2)
    @test_throws AssertionError QuanticsGrids.grid_origcoords(grid_2d, 3)
    @test_throws AssertionError QuanticsGrids.grid_origcoords(grid_custom, 4)

    # Test that returned range is iterable and has correct type
    coords = QuanticsGrids.grid_origcoords(grid_1d, 1)
    @test coords isa AbstractRange
    @test eltype(coords) <: AbstractFloat

    # Test consistency with grid spacing
    grid_spacing = QuanticsGrids.grid_step(grid_2d)
    coords_x_spacing = coords_x[2] - coords_x[1]
    coords_y_spacing = coords_y[2] - coords_y[1]
    @test coords_x_spacing ≈ grid_spacing[1]
    @test coords_y_spacing ≈ grid_spacing[2]

    # CHALLENGING TEST: Complex multi-dimensional grid with asymmetric resolutions
    complex_Rs = (5, 3, 7, 2)  # Different resolutions per dimension
    grid_complex = NewDiscretizedGrid(complex_Rs;
        lower_bound=(-10.0, 0.1, 50.0, -2.5),
        upper_bound=(15.0, 3.9, 100.0, 7.8),
        base=3)

    # Test each dimension has correct properties
    for d in 1:4
        coords_d = QuanticsGrids.grid_origcoords(grid_complex, d)
        expected_length = 3^complex_Rs[d]
        @test length(coords_d) == expected_length
        @test coords_d[1] ≈ QuanticsGrids.grid_min(grid_complex)[d]
        @test coords_d[end] ≈ QuanticsGrids.grid_max(grid_complex)[d]

        # Test uniform spacing within each dimension
        spacings = diff(collect(coords_d))
        @test all(spacing -> spacing ≈ spacings[1], spacings)
        @test spacings[1] ≈ QuanticsGrids.grid_step(grid_complex)[d]

        # Test consistency with grididx_to_origcoord
        test_indices = [1, expected_length ÷ 2, expected_length]
        for idx in test_indices
            grid_idx = ntuple(i -> i == d ? idx : 1, 4)
            expected_coord = grididx_to_origcoord(grid_complex, grid_idx)[d]
            @test coords_d[idx] ≈ expected_coord
        end
    end

    # Test with very small grid (edge case)
    grid_tiny = NewDiscretizedGrid{1}(1; lower_bound=(0.0,), upper_bound=(1.0,))
    coords_tiny = QuanticsGrids.grid_origcoords(grid_tiny, 1)
    @test length(coords_tiny) == 2
    @test coords_tiny[1] ≈ 0.0
    @test coords_tiny[2] ≈ 0.5  # For R=1, max coord should be 1 - 1/2 = 0.5
end
