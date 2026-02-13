@testitem "DiscretizedGrid constructors error paths" begin
    # upper_bound < lower_bound
    @test_throws ArgumentError DiscretizedGrid(10, 0.1, 0.01)

    # upper_bound < lower_bound (multiple) dimensions
    @test_throws ArgumentError DiscretizedGrid((10, 4), (0.1, 0.2), (0.01, 0.02))

    # upper_bound == lower_bound
    @test_throws ArgumentError DiscretizedGrid(10, 0.1, 0.1)
    @test_throws ArgumentError DiscretizedGrid((10, 4), (0.1, 0.2), (0.1, 0.02))

    # incompatible indextable/variablenames
    @test_throws ArgumentError DiscretizedGrid((:x, :y), [[(:x, 1)], [(:z, 1)]])

    # R = 0 and includeendpoint=true
    @test_throws ArgumentError DiscretizedGrid((0, 5); includeendpoint=true)
end

@testitem "DiscretizedGrid kwarg handling: _handle_kwargs_input" begin
    g = DiscretizedGrid((7, 8, 9); variablenames=(:x, :y, :z))
    quants_kwarg = origcoord_to_quantics(g; x=0.5, z=0.75, y=0.25)
    quants_plain = origcoord_to_quantics(g, (0.5, 0.25, 0.75))
    @test quants_kwarg == quants_plain
    @test_throws ArgumentError origcoord_to_quantics(g; x=0.5, z=0.75)  # missing y
    @test_throws ArgumentError origcoord_to_quantics(g; x=0.5, z=0.75, y=0.25, extra=1)  # extra keyword argument
    @test_throws ArgumentError origcoord_to_quantics(g; x=:abc, y=0.3, z=0.1)
    # @test_throws ArgumentError false
end

@testitem "DiscretizedGrid sitedim" begin
    g = DiscretizedGrid((7, 8, 9); variablenames=(:x, :y, :z))
    expected_sitedims = [8, 8, 8, 8, 8, 8, 8, 4, 2]
    siteinds = eachindex(QuanticsGrids.grid_indextable(g))
    for i in siteinds
        @test QuanticsGrids.sitedim(g, i) == expected_sitedims[i]
    end
    @test_throws DomainError QuanticsGrids.sitedim(g, last(siteinds) + 1)
end

@testitem "DiscretizedGrid constructors consistency" begin
    Rs = (2, 3, 4)
    variablenames = (:x, :y, :z)
    indextable = [[(:z, 1), (:y, 1), (:x, 1)], [(:z, 2), (:y, 2), (:x, 2)], [(:z, 3), (:y, 3)], [(:z, 4)]]
    g = DiscretizedGrid(Rs; variablenames)
    g´ = DiscretizedGrid(variablenames, Rs)
    g´´ = DiscretizedGrid(variablenames, indextable)
    @test g.discretegrid.lookup_table == g´.discretegrid.lookup_table == g´´.discretegrid.lookup_table
end

@testitem "DiscretizedGrid accepts untyped indextable container" begin
    R = 10
    N = 2^R
    lower_bound = (0.0, 0.0)
    upper_bound = (2.0, N - 1.0)

    indextable = []
    for l in 1:R
        push!(indextable, [(:w, l)])
        push!(indextable, [(:n, R - l + 1)])
    end

    g = DiscretizedGrid((:w, :n), indextable; lower_bound, upper_bound, includeendpoint=(false, true))
    @test QuanticsGrids.grid_variablenames(g) == (:w, :n)
    @test length(QuanticsGrids.grid_indextable(g)) == 2R
    @test QuanticsGrids.grid_max(g)[2] ≈ N - 1
end

@testitem "DiscretizedGrid dimension mismatch for bounds" begin
    @test_throws ArgumentError DiscretizedGrid{2}(3, (0.0,), (1.0, 2.0))
    @test_throws ArgumentError DiscretizedGrid{2}(3, (0.0, 1.0), (1.0,))
end

@testitem "DiscretizedGrid mixed-type bounds" begin
    R = 10
    N = 2^R
    g = DiscretizedGrid{2}((R, R), (0.0, 0), (2.0, N - 1); includeendpoint=(false, true))

    @test QuanticsGrids.lower_bound(g) == (0.0, 0.0)
    @test QuanticsGrids.grid_max(g)[2] ≈ N - 1
end

@testitem "DiscretizedGrid 0-dimensional show method" begin
    g = DiscretizedGrid(())
    @test try
        sprint(show, g)
        true
    catch e
        false
    end

    g = DiscretizedGrid((0,))
    @test try
        sprint(show, g)
        true
    catch e
        false
    end
end

@testitem "DiscretizedGrid show method" begin
    g = DiscretizedGrid((2, 3, 4))
    text = sprint(show, MIME"text/plain"(), g)
    @test occursin("Index table: [", text)
    @test occursin("1:(", text)
    @test !occursin("Symbol(\"", text)
    @test !occursin("\n│  ", text)
end
