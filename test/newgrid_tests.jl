using QuanticsGrids
using Test

@testset "NewDiscretizedGrid" begin
    @testset "constructor, square grid" begin
        grid = QuanticsGrids.NewDiscretizedGrid{2}(10)
        @test grid.Rs == (10, 10)
        @test only(grid.indextable[1]).variablename == Symbol(1)
        @test only(grid.indextable[1]).bitnumber == 1
        @test only(grid.indextable[2]).variablename == Symbol(2)
        @test only(grid.indextable[2]).bitnumber == 1
    end

    @testset "constructor, rectangular grid" begin
        grid = QuanticsGrids.NewDiscretizedGrid((3, 5))
        @test only(grid.indextable[6]).variablename == Symbol(2)
        @test only(grid.indextable[6]).bitnumber == 3
        @test only(grid.indextable[7]).variablename == Symbol(2)
        @test only(grid.indextable[7]).bitnumber == 4
        @test only(grid.indextable[8]).variablename == Symbol(2)
        @test only(grid.indextable[8]).bitnumber == 5
    end

    @testset "quantics_to_grididx, rectangular grid" begin
        grid = QuanticsGrids.NewDiscretizedGrid((3, 5))
        @test QuanticsGrids.quantics_to_grididx(grid, [1, 2, 1, 2, 1, 2, 1, 2]) == (1, 30)
    end

    @testset "grididx_to_quantics, rectangular grid" begin
        grid = QuanticsGrids.NewDiscretizedGrid((3, 5))
        @test QuanticsGrids.grididx_to_quantics(grid, (1, 30)) == [1, 2, 1, 2, 1, 2, 1, 2]
    end

    @testset "quantics_to_grididx ∘ grididx_to_quantics == identity" begin
        grid = QuanticsGrids.NewDiscretizedGrid((48, 31, 62); base=13)
        for _ in 1:100
            grididx = ntuple(d -> rand(1:grid.Rs[d]), ndims(grid))
            @test quantics_to_grididx(grid, grididx_to_quantics(grid, grididx)) == grididx
        end
    end

    @testset "grididx_to_quantics ∘ quantics_to_grididx == identity" begin
        grid = QuanticsGrids.NewDiscretizedGrid((48, 31, 62))
        for _ in 1:100
            quantics = rand(1:2, length(grid))
            @test grididx_to_quantics(grid, quantics_to_grididx(grid, quantics)) == quantics
        end
    end
end
