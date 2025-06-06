@testitem "NewDiscretizedGrid comprehensive functionality test" begin
    variablenames = (:x, :y, :z, :w)
    complex_tupletable = [
        [(:x, 1), (:y, 1), (:z, 1)],
        [(:w, 1)],
        [(:z, 2), (:w, 2), (:x, 3), (:y, 3)],
        [(:x, 2), (:y, 2)],
        [(:z, 3)],
        [(:w, 3), (:x, 4)],
        [(:y, 4), (:z, 4), (:w, 4)]
    ]

    lower_bound = (-π, 1e-12, 1e8, -1e6)
    upper_bound = (2 * π, 1e-10, 1e8 + 1000, 1e6)
    base = 3

    grid = QuanticsGrids.NewDiscretizedGrid(variablenames, complex_tupletable;
        lower_bound=lower_bound, upper_bound=upper_bound, base=base, includeendpoint=true)

    @test ndims(grid) == 4
    @test length(grid) == 7
    @test grid.base == 3
    @test grid.variablenames == variablenames

    expected_Rs = (4, 4, 4, 4)
    @test grid.Rs == expected_Rs

    expected_localdims = [27, 3, 81, 9, 3, 9, 27]
    @test QuanticsGrids.localdimensions(grid) == expected_localdims

    for (i, expected_dim) in enumerate(expected_localdims)
        @test QuanticsGrids.sitedim(grid, i) == expected_dim
    end

    lb = QuanticsGrids.lower_bound(grid)
    ub = QuanticsGrids.upper_bound(grid)
    gmin = QuanticsGrids.grid_min(grid)
    gmax = QuanticsGrids.grid_max(grid)
    gstep = QuanticsGrids.grid_step(grid)
    gorigin = QuanticsGrids.grid_origin(grid)

    @test lb == lower_bound
    @test gmin == lower_bound
    @test gorigin == lower_bound

    @test ub != upper_bound
    @test all(ub .> upper_bound)

    @test all(gmax .≈ upper_bound)

    expected_steps = (ub .- lb) ./ (base .^ grid.Rs)
    @test all(gstep .≈ expected_steps)

    corner_indices = [
        (1, 1, 1, 1),
        (base^4, base^4, base^4, base^4),
        (1, base^4, 1, base^4),
        (base^4, 1, base^4, 1),
        (base^2, base^2, base^2, base^2),
    ]

    for grididx in corner_indices
        quantics = grididx_to_quantics(grid, grididx)
        recovered_grididx = quantics_to_grididx(grid, quantics)
        @test recovered_grididx == grididx

        origcoord = grididx_to_origcoord(grid, grididx)
        recovered_grididx2 = QuanticsGrids.origcoord_to_grididx(grid, origcoord)
        @test recovered_grididx2 == grididx

        origcoord2 = quantics_to_origcoord(grid, quantics)
        recovered_quantics = origcoord_to_quantics(grid, origcoord2)
        @test recovered_quantics == quantics

        @test all(origcoord .≈ origcoord2)

        @test all(lb .<= origcoord .<= ub)
    end

    for _ in 1:50
        random_quantics = [rand(1:QuanticsGrids.sitedim(grid, site)) for site in 1:length(grid)]

        grididx = quantics_to_grididx(grid, random_quantics)
        origcoord = grididx_to_origcoord(grid, grididx)
        back_to_grididx = QuanticsGrids.origcoord_to_grididx(grid, origcoord)
        back_to_quantics = grididx_to_quantics(grid, back_to_grididx)
        back_to_origcoord = quantics_to_origcoord(grid, back_to_quantics)

        @test back_to_grididx == grididx
        @test back_to_quantics == random_quantics
        @test all(back_to_origcoord .≈ origcoord)

        @test all(1 .<= grididx .<= (base .^ grid.Rs))
        @test all(lb .<= origcoord .<= ub)
        @test length(back_to_quantics) == length(grid)
    end

    boundary_coords = [
        lb,
        ub,
        upper_bound,
        gmax,
    ]

    for coord in boundary_coords
        if all(lb .<= coord .<= ub)
            grididx = QuanticsGrids.origcoord_to_grididx(grid, coord)
            recovered_coord = grididx_to_origcoord(grid, grididx)
            re_recovered_grididx = QuanticsGrids.origcoord_to_grididx(grid, recovered_coord)
            @test re_recovered_grididx == grididx
        end
    end

    epsilon = 1e-10
    for i in 1:4
        coord_inside = ntuple(j -> j == i ? lb[j] + epsilon : (lb[j] + ub[j]) / 2, 4)
        if all(lb .<= coord_inside .<= ub)
            grididx = QuanticsGrids.origcoord_to_grididx(grid, coord_inside)
            @test all(1 .<= grididx .<= (base .^ grid.Rs))
        end
    end

    @test_throws AssertionError QuanticsGrids.origcoord_to_grididx(grid, (lb[1] - 1, lb[2], lb[3], lb[4]))
    @test_throws AssertionError QuanticsGrids.origcoord_to_grididx(grid, (ub[1] + 1, ub[2], ub[3], ub[4]))

    @test_throws AssertionError grididx_to_origcoord(grid, (0, 1, 1, 1))
    @test_throws AssertionError grididx_to_origcoord(grid, (base^4 + 1, 1, 1, 1))

    invalid_quantics = [28, 1, 1, 1, 1, 1, 1]
    @test_throws AssertionError quantics_to_grididx(grid, invalid_quantics)

    precision_test_coords = [
        (lb[1] + gstep[1] * 0.5, lb[2], lb[3], lb[4]),
        (lb[1], lb[2] + gstep[2] * 0.999999, lb[3], lb[4]),
        (ub[1] - gstep[1] * 1e-15, ub[2], ub[3], ub[4]),
    ]

    for coord in precision_test_coords
        if all(lb .<= coord .<= ub)
            grididx = QuanticsGrids.origcoord_to_grididx(grid, coord)
            recovered_coord = grididx_to_origcoord(grid, grididx)
            re_recovered_grididx = QuanticsGrids.origcoord_to_grididx(grid, recovered_coord)
            @test re_recovered_grididx == grididx
        end
    end

    test_function(x, y, z, w) = sin(x) * cos(y) * exp(-z / 1e8) * (w^2)
    quantics_func = QuanticsGrids.quanticsfunction(Float64, grid, test_function)

    for _ in 1:10
        random_quantics = [rand(1:QuanticsGrids.sitedim(grid, site)) for site in 1:length(grid)]

        result1 = quantics_func(random_quantics)

        origcoord = quantics_to_origcoord(grid, random_quantics)
        result2 = test_function(origcoord...)

        @test result1 ≈ result2
    end

    simple_grid = NewDiscretizedGrid((3, 4); base=5, lower_bound=(0.0, -1.0), upper_bound=(1.0, 1.0))
    @test simple_grid.Rs == (3, 4)
    @test simple_grid.base == 5
    @test QuanticsGrids.localdimensions(simple_grid) == [5, 5, 5, 5, 5, 5, 5]

    grid_1d = NewDiscretizedGrid{1}(8; base=2, lower_bound=(-5.0,), upper_bound=(10.0,))
    @test ndims(grid_1d) == 1
    @test grid_1d.Rs == (8,)

    for i in [1, 100, 256]
        coord = grididx_to_origcoord(grid_1d, i)
        back_idx = QuanticsGrids.origcoord_to_grididx(grid_1d, coord)
        @test back_idx == i
    end

    Rs_test = (3, 3)
    grid_fused = NewDiscretizedGrid(Rs_test; unfoldingscheme=:fused, base=2)
    grid_interleaved = NewDiscretizedGrid(Rs_test; unfoldingscheme=:interleaved, base=2)

    @test length(grid_fused) == 3
    @test length(grid_interleaved) == 6

    @test grid_fused.Rs == grid_interleaved.Rs

    test_grididx = (4, 6)
    coord_fused = grididx_to_origcoord(grid_fused, test_grididx)
    coord_interleaved = grididx_to_origcoord(grid_interleaved, test_grididx)
    @test all(coord_fused .≈ coord_interleaved)

    large_grid = NewDiscretizedGrid{2}(20; lower_bound=(0.0, 0.0), upper_bound=(1.0, 1.0))

    max_idx = 2^20
    extreme_indices = [1, max_idx ÷ 2, max_idx - 1, max_idx]
    for idx in extreme_indices
        grididx = (idx, idx)
        if all(1 .<= grididx .<= (2, 2) .^ large_grid.Rs)
            quantics = grididx_to_quantics(large_grid, grididx)
            recovered_grididx = quantics_to_grididx(large_grid, quantics)
            @test recovered_grididx == grididx
        end
    end

    Rs_comp = (4, 4)
    grid_no_endpoint = NewDiscretizedGrid(Rs_comp; includeendpoint=false)
    grid_with_endpoint = NewDiscretizedGrid(Rs_comp; includeendpoint=true)

    ub_no = QuanticsGrids.upper_bound(grid_no_endpoint)
    ub_with = QuanticsGrids.upper_bound(grid_with_endpoint)
    @test all(ub_with .> ub_no)

    gmax_no = QuanticsGrids.grid_max(grid_no_endpoint)
    gmax_with = QuanticsGrids.grid_max(grid_with_endpoint)
    @test all(gmax_with .> gmax_no)

    step_no = QuanticsGrids.grid_step(grid_no_endpoint)
    step_with = QuanticsGrids.grid_step(grid_with_endpoint)
    @test all(step_no .> 0) && all(step_with .> 0)

    @test isa(grid, NewDiscretizedGrid{4})
    @test all(isfinite.(QuanticsGrids.lower_bound(grid)))
    @test all(isfinite.(QuanticsGrids.upper_bound(grid)))
    @test all(QuanticsGrids.grid_step(grid) .> 0)

    @test length(grid.lookup_table) == 4
    @test all(length(grid.lookup_table[d]) == grid.Rs[d] for d in 1:4)

    base2_grid = NewDiscretizedGrid((5, 3); base=2, lower_bound=(0.0, -1.0), upper_bound=(3.0, 2.0))

    for _ in 1:20
        quantics_2 = [rand(1:QuanticsGrids.sitedim(base2_grid, site)) for site in 1:length(base2_grid)]
        grididx_2 = quantics_to_grididx(base2_grid, quantics_2)
        back_quantics_2 = grididx_to_quantics(base2_grid, grididx_2)
        @test back_quantics_2 == quantics_2

        origcoord_2 = grididx_to_origcoord(base2_grid, grididx_2)
        back_grididx_2 = QuanticsGrids.origcoord_to_grididx(base2_grid, origcoord_2)
        @test back_grididx_2 == grididx_2
    end

    simple_tupletable_base2 = [
        [(:x, 1)],
        [(:x, 2)],
        [(:y, 1)],
        [(:y, 2)]
    ]
    simple_grid_base2 = NewDiscretizedGrid((:x, :y), simple_tupletable_base2; base=2)

    test_quantics_base2 = [[1, 1, 1, 1], [2, 1, 1, 1], [1, 2, 1, 1], [2, 2, 2, 2]]
    for q in test_quantics_base2
        grididx = quantics_to_grididx(simple_grid_base2, q)
        back_q = grididx_to_quantics(simple_grid_base2, grididx)
        @test back_q == q
    end

    test_grid_2d = NewDiscretizedGrid((:x, :y), [[(:x, 1)], [(:y, 1)], [(:x, 2)], [(:y, 2)]];
        lower_bound=(-2.0, -3.0), upper_bound=(5.0, 7.0))

    grididx_kw = QuanticsGrids.origcoord_to_grididx(test_grid_2d; x=1.0, y=2.0)
    grididx_pos = QuanticsGrids.origcoord_to_grididx(test_grid_2d, (1.0, 2.0))
    @test grididx_kw == grididx_pos

    quantics_kw = origcoord_to_quantics(test_grid_2d; x=1.0, y=2.0)
    quantics_pos = origcoord_to_quantics(test_grid_2d, (1.0, 2.0))
    @test quantics_kw == quantics_pos

    origcoord_kw = grididx_to_origcoord(test_grid_2d; x=2, y=3)
    origcoord_pos = grididx_to_origcoord(test_grid_2d, (2, 3))
    @test all(origcoord_kw .≈ origcoord_pos)

    quantics_from_grididx_kw = grididx_to_quantics(test_grid_2d; x=2, y=3)
    quantics_from_grididx_pos = grididx_to_quantics(test_grid_2d, (2, 3))
    @test quantics_from_grididx_kw == quantics_from_grididx_pos

    test_grid_3d = NewDiscretizedGrid((:a, :b, :c), [[(:a, 1)], [(:b, 1)], [(:c, 1)]];
        lower_bound=(-1.0, 0.0, 10.0), upper_bound=(1.0, 2.0, 20.0))

    grididx_3d_kw = QuanticsGrids.origcoord_to_grididx(test_grid_3d; a=0.5, b=1.0, c=15.0)
    grididx_3d_pos = QuanticsGrids.origcoord_to_grididx(test_grid_3d, (0.5, 1.0, 15.0))
    @test grididx_3d_kw == grididx_3d_pos

    @test_throws AssertionError QuanticsGrids.origcoord_to_grididx(test_grid_2d; x=1.0)
    @test_throws AssertionError QuanticsGrids.origcoord_to_grididx(test_grid_2d; y=2.0)
    @test_throws AssertionError QuanticsGrids.origcoord_to_grididx(test_grid_2d; x=1.0, y=2.0, z=3.0)
    @test_throws AssertionError QuanticsGrids.origcoord_to_grididx(test_grid_2d; a=1.0, b=2.0)

    @test_throws AssertionError grididx_to_origcoord(test_grid_2d; x=1)
    @test_throws AssertionError grididx_to_origcoord(test_grid_2d; y=2)
    @test_throws AssertionError grididx_to_origcoord(test_grid_2d; x=1, y=2, z=3)
    @test_throws AssertionError grididx_to_origcoord(test_grid_2d; a=1, b=2)

    for _ in 1:10
        coord_vals = (rand() * 6 - 2, rand() * 9 - 3)
        grididx_tuple = QuanticsGrids.origcoord_to_grididx(test_grid_2d, coord_vals)
        local grididx_kw = QuanticsGrids.origcoord_to_grididx(test_grid_2d; x=coord_vals[1], y=coord_vals[2])
        @test grididx_tuple == grididx_kw

        quantics_tuple = origcoord_to_quantics(test_grid_2d, coord_vals)
        local quantics_kw = origcoord_to_quantics(test_grid_2d; x=coord_vals[1], y=coord_vals[2])
        @test quantics_tuple == quantics_kw

        origcoord_from_grid_tuple = grididx_to_origcoord(test_grid_2d, grididx_tuple)
        origcoord_from_grid_kw = grididx_to_origcoord(test_grid_2d; x=grididx_tuple[1], y=grididx_tuple[2])
        @test all(origcoord_from_grid_tuple .≈ origcoord_from_grid_kw)

        quantics_from_grid_tuple = grididx_to_quantics(test_grid_2d, grididx_tuple)
        quantics_from_grid_kw = grididx_to_quantics(test_grid_2d; x=grididx_tuple[1], y=grididx_tuple[2])
        @test quantics_from_grid_tuple == quantics_from_grid_kw
    end
end
