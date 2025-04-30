
@testitem "grid.jl" begin
    using Test
    import QuanticsGrids as QD


    @testset "quanticsfunction" begin
        R = 8
        grid = QD.DiscretizedGrid{1}(R, 0.0, 1.0)
        fx(x) = exp(-x)
        fq = QD.quanticsfunction(Float64, grid, fx)

        @test fq(ones(Int, R)) == fx(0.0)
        @test fq(fill(2, R)) ≈ fx(1.0 - 1 / 2^R)
    end

    @testset "1D grid (large R)" begin
        R = 62
        grid = QD.DiscretizedGrid{1}(R, 0.0, 1.0)
        @test QD.grididx_to_quantics(grid, 2^R) == fill(2, R)
    end

    @testset "2D grid (large R)" for base in [2]
        R = 62
        d = 2
        grid = QD.DiscretizedGrid{d}(R, 0.0, 1.0; base=base)
        @test QD.grididx_to_quantics(grid, ntuple(i -> base^R, d)) == fill(base^d, R)
    end

    @testset "1D grid (too large R)" begin
        R = 64
        @test_throws ErrorException QD.DiscretizedGrid{1}(R, 0.0, 1.0)
    end

    @testset "grid representation conversion" for R in [10]
        reprs = [:grididx, :quantics, :origcoord]

        testset = [
            (QD.DiscretizedGrid{1}(R, 0.0, 1.0), 2),
            (QD.DiscretizedGrid{2}(R, (0.0, 0.0), (1.0, 1.0)), (2, 3)),
            (QD.InherentDiscreteGrid{1}(R), 2),
            (QD.InherentDiscreteGrid{2}(R), (2, 3)),
        ]
        for (grid, ini_grididx) in testset

            data = Dict{Symbol,Any}()
            transforms = Dict(
                (:grididx, :quantics) => QD.grididx_to_quantics,
                (:grididx, :origcoord) => QD.grididx_to_origcoord,
                (:quantics, :grididx) => QD.quantics_to_grididx,
                (:quantics, :origcoord) => QD.quantics_to_origcoord,
                (:origcoord, :grididx) => QD.origcoord_to_grididx,
                (:origcoord, :quantics) => QD.origcoord_to_quantics,
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

    @testset "InherentDiscreteGrid" for unfoldingscheme in [:interleaved, :fused],
        step in [(1, 1, 1), (1, 1, 2)],
        origin in [(1, 1, 1), (1, 1, 2)]

        m = QuanticsGrids.InherentDiscreteGrid{3}(5, origin; unfoldingscheme, step=step)
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

    @testset "DiscretizedGrid" for unfoldingscheme in [:interleaved, :fused]
        @testset "1D" begin
            R = 5
            grid_min = 0.1
            grid_max = 2.0
            dx = (grid_max - grid_min) / 2^R
            g = QuanticsGrids.DiscretizedGrid{1}(R, grid_min, grid_max; unfoldingscheme)
            @test QuanticsGrids.localdimensions(g) == fill(2, R)

            @test @inferred(
                QuanticsGrids.origcoord_to_grididx(g, 0.999999 * dx + grid_min)
            ) == 2
            @test QuanticsGrids.origcoord_to_grididx(g, 1.999999 * dx + grid_min) == 3
            @test QuanticsGrids.origcoord_to_grididx(g, grid_max - 1e-9 * dx - dx) == 2^R
            @test QuanticsGrids.grid_min(g) == 0.1
            @test QuanticsGrids.grid_max(g) == 2.0
            @test QuanticsGrids.grid_step(g) == 0.059375
        end

        @testset "1D (includeendpoint)" for unfoldingscheme in [:interleaved, :fused]
            R = 5
            grid_min = 0.0
            grid_max = 1.0
            dx = (grid_max - grid_min) / (2^R - 1)
            g = QuanticsGrids.DiscretizedGrid{1}(
                R,
                grid_min,
                grid_max;
                unfoldingscheme,
                includeendpoint=true,
            )
            @test QuanticsGrids.localdimensions(g) == fill(2, R)

            @test @inferred(QuanticsGrids.origcoord_to_grididx(g, grid_min)) == 1
            @test @inferred(QuanticsGrids.origcoord_to_grididx(g, grid_max)) == 2^R
            @test only(QuanticsGrids.grid_step(g)) == dx
            @test only(QuanticsGrids.quantics_to_origcoord(g, fill(2, R))) == grid_max
        end

        @testset "2D" for unfoldingscheme in [:interleaved, :fused]
            R = 5
            d = 2
            grid_min = (0.1, 0.1)
            grid_max = (2.0, 2.0)
            dx = (grid_max .- grid_min) ./ 2^R
            g = QuanticsGrids.DiscretizedGrid{d}(R, grid_min, grid_max; unfoldingscheme)

            if unfoldingscheme === :interleaved
                @test QuanticsGrids.localdimensions(g) == fill(2, d * R)
            else
                @test QuanticsGrids.localdimensions(g) == fill(2^d, R)
            end

            @test QuanticsGrids.grid_min(g) == (0.1, 0.1)
            @test QuanticsGrids.grid_step(g) == dx == (0.059375, 0.059375)
            @test QuanticsGrids.grid_max(g) == (2.0, 2.0)

            cs = [
                0.999999 .* dx .+ grid_min,
                1.999999 .* dx .+ grid_min,
                grid_max .- 1e-9 .* dx .- dx,
            ]
            refs = [2, 3, 2^R]

            for (c, ref) in zip(cs, refs)
                @inferred(QuanticsGrids.origcoord_to_grididx(g, c))
                @test all(QuanticsGrids.origcoord_to_grididx(g, c) .== ref)
            end

            @test_throws "Bound Error:" QuanticsGrids.origcoord_to_grididx(g, (0.0, 0.0))
            @test_throws "Bound Error:" QuanticsGrids.origcoord_to_grididx(g, (0.0, 1.1))
            @test_throws "Bound Error:" QuanticsGrids.origcoord_to_grididx(g, (1.1, 0.0))
            @test_throws "Bound Error:" QuanticsGrids.origcoord_to_grididx(g, (3.0, 1.1))
            @test_throws "Bound Error:" QuanticsGrids.origcoord_to_grididx(g, (1.1, 3.0))
            @test_throws "Bound Error:" QuanticsGrids.origcoord_to_grididx(g, (3.0, 3.0))
        end

        @testset "2D (includeendpoint)" for unfoldingscheme in [:interleaved, :fused]
            R = 5
            d = 2
            grid_min = (0.1, 0.1)
            grid_max = (2.0, 2.0)
            dx = (grid_max .- grid_min) ./ (2^R - 1)
            g = QuanticsGrids.DiscretizedGrid{d}(
                R,
                grid_min,
                grid_max;
                unfoldingscheme,
                includeendpoint=true,
            )
            @test g.includeendpoint
            @test QuanticsGrids.grid_step(g) == dx

            if unfoldingscheme === :interleaved
                @test QuanticsGrids.localdimensions(g) == fill(2, d * R)
            else
                @test QuanticsGrids.localdimensions(g) == fill(2^d, R)
            end
        end
    end

    @testset "dimension inference constructors" begin
        g1 = QD.DiscretizedGrid(4, 0.0, 1.0)
        @test g1 isa QD.DiscretizedGrid{1}
        @test g1.grid_min == (0.0,)
        @test g1.grid_max == (1.0,)

        g2 = QD.DiscretizedGrid(4, (0.0, 0.0), (1.0, 1.0))
        @test g2 isa QD.DiscretizedGrid{2}
        @test g2.grid_min == (0.0, 0.0)
        @test g2.grid_max == (1.0, 1.0)

        g6 = QD.DiscretizedGrid(6, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0), (1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
        @test g6 isa QD.DiscretizedGrid{6}
        @test g6.grid_min == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        @test g6.grid_max == (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

        g3 = QD.InherentDiscreteGrid(4, 1)
        @test g3 isa QD.InherentDiscreteGrid{1}
        @test g3.origin == (1,)

        g4 = QD.InherentDiscreteGrid(4, (1, 2, 3); step=(1, 2, 1))
        @test g4 isa QD.InherentDiscreteGrid{3}
        @test g4.origin == (1, 2, 3)
        @test g4.step == (1, 2, 1)
    end
end
