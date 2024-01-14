
@testitem "grid.jl" begin
    using Test
    import QuanticsGrids as QD

    @testset "grid representation conversion" begin
        reprs = [:grididx, :quantics, :origcoord]

        R = 10
        testset = [
            (QD.DiscretizedGrid{1}(R, (0.0,), (1.0,)), (2,)),
            (QD.DiscretizedGrid{2}(R, (0.0,0.0), (1.0,1.0)), (2,3)),
            (QD.InherentDiscreteGrid{1}(R), (2,)),
            (QD.InherentDiscreteGrid{2}(R), (2,3)),
        ]
        for (grid, initial) in testset

            data = Dict{Symbol,Any}()
            transforms = Dict(
                (:grididx, :quantics) => QD.grididx_to_quantics,
                (:grididx, :origcoord) => QD.grididx_to_origcoord,
                (:quantics, :grididx) => QD.quantics_to_grididx,
                (:quantics, :origcoord) => QD.quantics_to_origcoord,
                (:origcoord, :grididx) => QD.origcoord_to_grididx,
                (:origcoord, :quantics) => QD.origcoord_to_quantics,
            )

            data[:grididx] = initial
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

    @testset "InherentDiscreteGrid" for unfoldingscheme in instances(QuanticsGrids.UnfoldingSchemes.UnfoldingScheme)
        m = QuanticsGrids.InherentDiscreteGrid{3}(5; unfoldingscheme)
        @test QuanticsGrids.grid_min(m) == (1, 1, 1)
        @test QuanticsGrids.grid_step(m) == (1, 1, 1)
        for idx in [(1, 1, 1), (1, 1, 2), (1, 25, 1), (14, 1, 1), (25, 25, 25)]
            c = QuanticsGrids.grididx_to_origcoord(m, idx)
            @test QuanticsGrids.origcoord_to_grididx(m, c) == idx

            q = QuanticsGrids.grididx_to_quantics(m, idx)
            if unfoldingscheme == QuanticsGrids.UnfoldingSchemes.fused
                @test length(q) == 5
            else
                @test length(q) == 3 * 5
            end
            @test all((1 .<= q) .&& (q .<= 2^3))
            @test QuanticsGrids.quantics_to_origcoord(m, q) == c
        end
    end

    @testset "DiscretizedGrid" for unfoldingscheme in instances(QuanticsGrids.UnfoldingSchemes.UnfoldingScheme)
        @testset "1D" begin
            R = 5
            grid_min = 0.1
            grid_max = 2.0
            dx = (grid_max - grid_min) / 2^R
            g = QuanticsGrids.DiscretizedGrid{1}(R, (grid_min,), (grid_max,); unfoldingscheme)
            @test @inferred(
                QuanticsGrids.origcoord_to_grididx(g, 0.999999 * dx + grid_min)
            ) == 1
            @test QuanticsGrids.origcoord_to_grididx(g, 1.999999 * dx + grid_min) == 2
            @test QuanticsGrids.origcoord_to_grididx(g, grid_max - 1e-9 * dx) == 2^R
            @test QuanticsGrids.grid_min(g) == (0.1,)
            @test QuanticsGrids.grid_max(g) == (2.0,)
            @test QuanticsGrids.grid_step(g) == (0.059375,)
        end

        @testset "2D" for unfoldingscheme in instances(QuanticsGrids.UnfoldingSchemes.UnfoldingScheme)
            R = 5
            d = 2
            grid_min = (0.1, 0.1)
            grid_max = (2.0, 2.0)
            dx = (grid_max .- grid_min) ./ 2^R
            g = QuanticsGrids.DiscretizedGrid{d}(R, grid_min, grid_max; unfoldingscheme)

            @test QuanticsGrids.grid_min(g) == (0.1, 0.1)
            @test QuanticsGrids.grid_step(g) == dx == (0.059375, 0.059375)
            @test QuanticsGrids.grid_max(g) == (2.0, 2.0)

            cs = [
                0.999999 .* dx .+ grid_min,
                1.999999 .* dx .+ grid_min,
                grid_max .- 1e-9 .* dx,
            ]
            refs = [1, 2, 2^R]

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
    end
end
