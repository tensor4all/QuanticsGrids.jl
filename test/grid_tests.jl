@testitem "grid.jl" begin
    using Test
    import QuanticsGrids

    @testset "quanticsfunction" begin
        R = 8
        grid = QuanticsGrids.DiscretizedGrid{1}(R, 0.0, 1.0)
        fx(x) = exp(-x)
        fq = QuanticsGrids.quanticsfunction(Float64, grid, fx)

        @test fq(ones(Int, R)) == fx(0.0)
        @test fq(fill(2, R)) ≈ fx(1.0 - 1 / 2^R)
    end

    @testset "1D grid (large R)" begin
        R = 62
        grid = QuanticsGrids.DiscretizedGrid{1}(R, 0.0, 1.0)
        @test QuanticsGrids.grididx_to_quantics(grid, 2^R) == fill(2, R)
    end

    @testset "2D grid (large R)" for base in [2]
        R = 62
        d = 2
        grid = QuanticsGrids.DiscretizedGrid{d}(R, 0.0, 1.0; base=base)
        @test QuanticsGrids.grididx_to_quantics(grid, ntuple(i -> base^R, d)) ==
              fill(base^d, R)
    end

    @testset "1D grid (too large R)" begin
        R = 64
        @test_throws ArgumentError QuanticsGrids.DiscretizedGrid{1}(R, 0.0, 1.0)
    end

    @testset "grid representation conversion" for R in [10]
        reprs = [:grididx, :quantics, :origcoord]

        testset = [
            (QuanticsGrids.DiscretizedGrid{1}(R, 0.0, 1.0), 2),
            (QuanticsGrids.DiscretizedGrid{2}(R, (0.0, 0.0), (1.0, 1.0)), (2, 3)),
            (QuanticsGrids.InherentDiscreteGrid{1}(R), 2),
            (QuanticsGrids.InherentDiscreteGrid{2}(R), (2, 3)),
        ]
        for (grid, ini_grididx) in testset

            data = Dict{Symbol,Any}()
            transforms = Dict(
                (:grididx, :quantics) => QuanticsGrids.grididx_to_quantics,
                (:grididx, :origcoord) => QuanticsGrids.grididx_to_origcoord,
                (:quantics, :grididx) => QuanticsGrids.quantics_to_grididx,
                (:quantics, :origcoord) => QuanticsGrids.quantics_to_origcoord,
                (:origcoord, :grididx) => QuanticsGrids.origcoord_to_grididx,
                (:origcoord, :quantics) => QuanticsGrids.origcoord_to_quantics,
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

    @testset "InherentDiscreteGrid" for
        unfoldingscheme in [:interleaved, :fused, :serial, :meander],
        step in [(1, 1, 1), (1, 1, 2)],
        origin in [(1, 1, 1), (1, 1, 2)]

        m = QuanticsGrids.InherentDiscreteGrid{3}(5, origin; unfoldingscheme, step=step)
        @test QuanticsGrids.grid_min(m) == origin
        @test QuanticsGrids.grid_step(m) == step

        if unfoldingscheme === :fused
            @test QuanticsGrids.localdimensions(m) == fill(2^3, 5)
        else
            @test QuanticsGrids.localdimensions(m) == fill(2, 3 * 5)
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

    @testset "DiscretizedGrid" for unfoldingscheme in [:interleaved, :fused]
        @testset "1D" begin
            unfoldingscheme = :interleaved

            R = 5
            a = 0.1
            b = 2.0
            dx = (b - a) / 2^R
            g = QuanticsGrids.DiscretizedGrid{1}(R, a, b; unfoldingscheme)
            @test QuanticsGrids.localdimensions(g) == fill(2, R)

            @test @inferred(QuanticsGrids.origcoord_to_grididx(g, 0.999999 * dx + a)) == 2
            @test QuanticsGrids.origcoord_to_grididx(g, 1.999999 * dx + a) == 3
            @test QuanticsGrids.origcoord_to_grididx(
                g,
                QuanticsGrids.grid_max(g) + 1e-9 * dx,
            ) == 2^R

            @test QuanticsGrids.lower_bound(g) == 0.1
            @test QuanticsGrids.upper_bound(g) == 2.0
            @test QuanticsGrids.grid_min(g) == a
            @test QuanticsGrids.grid_max(g) == b - dx

            @test QuanticsGrids.grid_step(g) == 0.059375
        end

        @testset "1D (includeendpoint)" for unfoldingscheme in [:interleaved, :fused]
            R = 5
            a = 0.0
            b = 1.0
            dx = (b - a) / (2^R - 1)
            g = QuanticsGrids.DiscretizedGrid{1}(
                R,
                a,
                b;
                unfoldingscheme,
                includeendpoint=true,
            )
            @test QuanticsGrids.localdimensions(g) == fill(2, R)

            @test QuanticsGrids.lower_bound(g) == 0.0
            @test QuanticsGrids.upper_bound(g) ≈ 1.0 + dx
            @test QuanticsGrids.grid_min(g) == a
            @test QuanticsGrids.grid_max(g) == b

            @test @inferred(QuanticsGrids.origcoord_to_grididx(g, a)) == 1
            @test @inferred(QuanticsGrids.origcoord_to_grididx(g, b)) == 2^R
            @test only(QuanticsGrids.grid_step(g)) == dx
            @test only(QuanticsGrids.quantics_to_origcoord(g, fill(2, R))) == b
        end

        @testset "2D" for unfoldingscheme in [:interleaved, :fused, :serial, :meander]
            R = 5
            d = 2
            a = (0.1, 0.1)
            b = (2.0, 2.0)
            dx = (b .- a) ./ 2^R
            g = QuanticsGrids.DiscretizedGrid{d}(R, a, b; unfoldingscheme)

            if unfoldingscheme === :fused
                @test QuanticsGrids.localdimensions(g) == fill(2^d, R)
            else
                @test QuanticsGrids.localdimensions(g) == fill(2, d * R)
            end

            @test QuanticsGrids.lower_bound(g) == (0.1, 0.1)
            @test QuanticsGrids.grid_step(g) == dx == (0.059375, 0.059375)
            @test QuanticsGrids.upper_bound(g) == (2.0, 2.0)

            @test QuanticsGrids.grid_min(g) == a
            @test QuanticsGrids.grid_max(g) == b .- dx

            cs = [0.999999 .* dx .+ a, 1.999999 .* dx .+ a, b .- dx .- 1e-9 .* dx]
            refs = [2, 3, 2^R]

            for (c, ref) in zip(cs, refs)
                @inferred(QuanticsGrids.origcoord_to_grididx(g, c))
                @test all(QuanticsGrids.origcoord_to_grididx(g, c) .== ref)
            end

            @test_throws DomainError QuanticsGrids.origcoord_to_grididx(g, (0.0, 0.0))
            @test_throws DomainError QuanticsGrids.origcoord_to_grididx(g, (0.0, 1.1))
            @test_throws DomainError QuanticsGrids.origcoord_to_grididx(g, (1.1, 0.0))
            @test_throws DomainError QuanticsGrids.origcoord_to_grididx(g, (3.0, 1.1))
            @test_throws DomainError QuanticsGrids.origcoord_to_grididx(g, (1.1, 3.0))
            @test_throws DomainError QuanticsGrids.origcoord_to_grididx(g, (3.0, 3.0))
        end

        @testset "2D (includeendpoint)" for unfoldingscheme in [:interleaved, :fused, :serial, :meander]
            R = 5
            d = 2
            a = (0.1, 0.1)
            b = (2.0, 2.0)
            dx = (b .- a) ./ (2^R - 1)
            g = QuanticsGrids.DiscretizedGrid{d}(
                R,
                a,
                b;
                unfoldingscheme,
                includeendpoint=true,
            )
            @test all(QuanticsGrids.grid_step(g) .≈ dx)

            if unfoldingscheme === :fused
                @test QuanticsGrids.localdimensions(g) == fill(2^d, R)
            else
                @test QuanticsGrids.localdimensions(g) == fill(2, d * R)
            end
        end
    end

    @testset "dimension inference constructors" begin
        g1 = QuanticsGrids.DiscretizedGrid(4, 0.0, 1.0)
        g1_ref = QuanticsGrids.DiscretizedGrid{1}(4, 0.0, 1.0)
        @test g1 isa QuanticsGrids.DiscretizedGrid{1}
        @test QuanticsGrids.grid_min(g1) == QuanticsGrids.grid_min(g1_ref)
        @test QuanticsGrids.grid_max(g1) == QuanticsGrids.grid_max(g1_ref)
        @test QuanticsGrids.upper_bound(g1) == QuanticsGrids.upper_bound(g1_ref)

        g2 = QuanticsGrids.DiscretizedGrid(4, (0.0, 0.0), (1.0, 1.0))
        g2_ref = QuanticsGrids.DiscretizedGrid{2}(4, (0.0, 0.0), (1.0, 1.0))
        @test g2 isa QuanticsGrids.DiscretizedGrid{2}
        @test QuanticsGrids.grid_min(g2) == QuanticsGrids.grid_min(g2_ref)
        @test QuanticsGrids.grid_max(g2) == QuanticsGrids.grid_max(g2_ref)
        @test QuanticsGrids.upper_bound(g2) == QuanticsGrids.upper_bound(g2_ref)

        g6 = QuanticsGrids.DiscretizedGrid(
            6,
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
        )
        g6_ref = QuanticsGrids.DiscretizedGrid{6}(
            6,
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
        )
        @test g6 isa QuanticsGrids.DiscretizedGrid{6}
        @test QuanticsGrids.grid_min(g6) == QuanticsGrids.grid_min(g6_ref)
        @test QuanticsGrids.grid_max(g6) == QuanticsGrids.grid_max(g6_ref)
        @test QuanticsGrids.upper_bound(g6) == QuanticsGrids.upper_bound(g6_ref)

        g3 = QuanticsGrids.InherentDiscreteGrid(4, 1)
        @test g3 isa QuanticsGrids.InherentDiscreteGrid{1}
        @test QuanticsGrids.grid_origin(g3) == 1

        g4 = QuanticsGrids.InherentDiscreteGrid(4, (1, 2, 3); step=(1, 2, 1))
        @test g4 isa QuanticsGrids.InherentDiscreteGrid{3}
        @test QuanticsGrids.grid_origin(g4) == (1, 2, 3)
        @test QuanticsGrids.grid_step(g4) == (1, 2, 1)
    end
end

@testitem "unfolding schemes :serial and :meander" begin
    using Test
    import QuanticsGrids

    schemes = [:fused, :interleaved, :serial, :meander]

    # 1D grids return a scalar index instead of a 1-tuple
    expected_index(d, grididx) = d == 1 ? grididx[1] : grididx

    @testset "roundtrip ($d dims, base = $base, R = $R, $unfoldingscheme)" for
        d in 1:3, base in (2, 3), R in 1:3, unfoldingscheme in schemes

        grid = QuanticsGrids.InherentDiscreteGrid{d}(
            ntuple(_ -> R, d); base, unfoldingscheme
        )
        for cartesian in CartesianIndices(ntuple(_ -> Base.OneTo(base^R), d))
            grididx = Tuple(cartesian)
            quantics = QuanticsGrids.grididx_to_quantics(grid, grididx)
            @test QuanticsGrids.quantics_to_grididx(grid, quantics) == expected_index(d, grididx)
        end
    end

    @testset "DiscretizedGrid roundtrip ($d dims, base = $base, R = $R, $unfoldingscheme)" for
        d in 1:2, base in (2, 3), R in 1:2, unfoldingscheme in schemes

        grid = QuanticsGrids.DiscretizedGrid{d}(R; base, unfoldingscheme)
        for cartesian in CartesianIndices(ntuple(_ -> Base.OneTo(base^R), d))
            grididx = Tuple(cartesian)
            quantics = QuanticsGrids.grididx_to_quantics(grid, grididx)
            @test QuanticsGrids.quantics_to_grididx(grid, quantics) == expected_index(d, grididx)
        end
    end

    # Pinned digit sequences: the base-2 digits of grididx - 1 are written most
    # significant first, i.e. X = 2 is x = (1, 2) and Y = 4 is y = (2, 2).
    @testset "digit ordering (d = 2, R = 1, base = 2)" begin
        grid_for(scheme) = QuanticsGrids.InherentDiscreteGrid{2}((1, 1); unfoldingscheme=scheme)
        for scheme in [:interleaved, :serial, :meander]
            @test QuanticsGrids.grididx_to_quantics(grid_for(scheme), (1, 1)) == [1, 1]
            @test QuanticsGrids.grididx_to_quantics(grid_for(scheme), (2, 1)) == [2, 1]
            @test QuanticsGrids.grididx_to_quantics(grid_for(scheme), (1, 2)) == [1, 2]
            @test QuanticsGrids.grididx_to_quantics(grid_for(scheme), (2, 2)) == [2, 2]
        end
        @test QuanticsGrids.grididx_to_quantics(grid_for(:fused), (1, 1)) == [1]
        @test QuanticsGrids.grididx_to_quantics(grid_for(:fused), (2, 1)) == [2]
        @test QuanticsGrids.grididx_to_quantics(grid_for(:fused), (1, 2)) == [3]
        @test QuanticsGrids.grididx_to_quantics(grid_for(:fused), (2, 2)) == [4]
    end

    # With more than one digit per variable the schemes differ: :serial is
    # x_1 x_2 y_1 y_2, :meander reverses the odd-numbered variables (x_2 x_1 y_1 y_2),
    # :interleaved is x_1 y_1 x_2 y_2 and :fused fuses (x_i, y_i) per digit level.
    @testset "digit ordering (d = 2, R = 2, base = 2)" begin
        grid_for(scheme) = QuanticsGrids.InherentDiscreteGrid{2}((2, 2); unfoldingscheme=scheme)

        @test QuanticsGrids.grididx_to_quantics(grid_for(:serial), (2, 1)) == [1, 2, 1, 1]
        @test QuanticsGrids.grididx_to_quantics(grid_for(:serial), (3, 4)) == [2, 1, 2, 2]
        @test QuanticsGrids.grididx_to_quantics(grid_for(:meander), (2, 1)) == [2, 1, 1, 1]
        @test QuanticsGrids.grididx_to_quantics(grid_for(:meander), (3, 4)) == [1, 2, 2, 2]
        @test QuanticsGrids.grididx_to_quantics(grid_for(:interleaved), (2, 1)) == [1, 1, 2, 1]
        @test QuanticsGrids.grididx_to_quantics(grid_for(:interleaved), (3, 4)) == [2, 2, 1, 2]
        @test QuanticsGrids.grididx_to_quantics(grid_for(:fused), (2, 1)) == [1, 2]
        @test QuanticsGrids.grididx_to_quantics(grid_for(:fused), (3, 4)) == [4, 3]

        @test QuanticsGrids.localdimensions(grid_for(:serial)) == [2, 2, 2, 2]
        @test QuanticsGrids.localdimensions(grid_for(:meander)) == [2, 2, 2, 2]
    end

    @testset ":serial is the :grouped ordering" for d in 1:3
        grid_for(scheme) = QuanticsGrids.InherentDiscreteGrid{d}(
            ntuple(_ -> 3, d); unfoldingscheme=scheme
        )
        @test QuanticsGrids.grid_indextable(grid_for(:serial)) ==
              QuanticsGrids.grid_indextable(grid_for(:grouped))
        if d > 1
            @test QuanticsGrids.grid_indextable(grid_for(:meander)) !=
                  QuanticsGrids.grid_indextable(grid_for(:serial))
        end
    end

    @testset "invalid unfoldingscheme is rejected" begin
        @test_throws ArgumentError QuanticsGrids.InherentDiscreteGrid{2}((3, 3); unfoldingscheme=:invalid)
        @test_throws ArgumentError QuanticsGrids.DiscretizedGrid{2}((3, 3); unfoldingscheme=:invalid)
    end
end
