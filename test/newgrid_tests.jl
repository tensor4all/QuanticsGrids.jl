@testitem "constructor, square grid" begin
    grid = QuanticsGrids.NewDiscretizedGrid{2}(10)
    @test grid.Rs == (10, 10)
    @test only(grid.indextable[1]).variablename == Symbol(1)
    @test only(grid.indextable[1]).bitnumber == 1
    @test only(grid.indextable[2]).variablename == Symbol(2)
    @test only(grid.indextable[2]).bitnumber == 1
end

@testitem "constructor, rectangular grid" begin
    grid = QuanticsGrids.NewDiscretizedGrid((3, 5))
    @test only(grid.indextable[6]).variablename == Symbol(2)
    @test only(grid.indextable[6]).bitnumber == 3
    @test only(grid.indextable[7]).variablename == Symbol(2)
    @test only(grid.indextable[7]).bitnumber == 4
    @test only(grid.indextable[8]).variablename == Symbol(2)
    @test only(grid.indextable[8]).bitnumber == 5
end

@testitem "quantics_to_grididx, rectangular grid" begin
    grid = QuanticsGrids.NewDiscretizedGrid((3, 5))
    @test QuanticsGrids.quantics_to_grididx(grid, [1, 2, 1, 2, 1, 2, 1, 2]) == (1, 30)
end

@testitem "grididx_to_quantics, rectangular grid" begin
    grid = QuanticsGrids.NewDiscretizedGrid((3, 5))
    @test QuanticsGrids.grididx_to_quantics(grid, (1, 30)) == [1, 2, 1, 2, 1, 2, 1, 2]
end

@testitem "quantics_to_grididx ∘ grididx_to_quantics == identity" begin
    grid = QuanticsGrids.NewDiscretizedGrid((5, 3, 17); base=13)
    for _ in 1:100
        grididx = ntuple(d -> rand(1:grid.Rs[d]), ndims(grid))
        @test QuanticsGrids.quantics_to_grididx(grid, QuanticsGrids.grididx_to_quantics(grid, grididx)) == grididx
    end
end

@testitem "grididx_to_quantics ∘ quantics_to_grididx == identity" begin
    grid = QuanticsGrids.NewDiscretizedGrid((48, 31, 62))
    for _ in 1:100
        quantics = rand(1:2, length(grid))
        @test QuanticsGrids.grididx_to_quantics(grid, QuanticsGrids.quantics_to_grididx(grid, quantics)) == quantics
    end
end

@testitem "grididx_to_quantics ∘ quantics_to_grididx == identity, base != 2" begin
    base = 7
    grid = QuanticsGrids.NewDiscretizedGrid((22, 9, 14); base)
    for _ in 1:100
        quantics = rand(1:base, length(grid))
        @test QuanticsGrids.grididx_to_quantics(grid, QuanticsGrids.quantics_to_grididx(grid, quantics)) == quantics
    end
end

@testitem "ctor from indextable" begin
    grid = QuanticsGrids.NewDiscretizedGrid((:a, :b, :c), [[(:a, 4)], [(:a, 3)], [(:a, 2)], [(:a, 1)], [(:b, 1)], [(:b, 2)], [(:b, 3)], [(:c, 1)], [(:c, 2)], [(:c, 3)]])
    @test grid.Rs == (4, 3, 3)
    @test grid.lower_bound == (0.0, 0.0, 0.0)
    @test grid.upper_bound == (1.0, 1.0, 1.0)
    @test grid.variablenames == (:a, :b, :c)
end

@testitem "ctor from indextable, quantics <-> grididx" begin
    grid = QuanticsGrids.NewDiscretizedGrid((:a, :b, :c), [[(:a, 4)], [(:a, 3)], [(:a, 2)], [(:a, 1)], [(:b, 1)], [(:b, 2)], [(:b, 3)], [(:c, 1)], [(:c, 2)], [(:c, 3)]])
    @test QuanticsGrids.quantics_to_grididx(grid, [1, 2, 1, 2, 1, 2, 1, 2, 1, 2]) == (11, 3, 6)
    @test QuanticsGrids.grididx_to_quantics(grid, (11, 3, 6)) == [1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
end

@testitem "ctor from indextable, quantics <-> grididx, fused indices" begin
    grid = QuanticsGrids.NewDiscretizedGrid((:a, :b, :c, :d), [[(:a, 4)], [(:a, 3)], [(:a, 2)], [(:a, 1)], [(:b, 1), (:d, 1)], [(:b, 2)], [(:b, 3)], [(:c, 1), (:d, 2)], [(:c, 2)], [(:c, 3)]])
    @test QuanticsGrids.quantics_to_grididx(grid, [1, 2, 1, 2, 2, 2, 1, 4, 1, 2]) == (11, 3, 6, 4)
    @test QuanticsGrids.grididx_to_quantics(grid, (11, 3, 6, 4)) == [1, 2, 1, 2, 2, 2, 1, 4, 1, 2]
end

@testitem "challenging tests - extreme edge cases" begin
    # Test with minimum valid grididx (all 1s)
    grid = QuanticsGrids.NewDiscretizedGrid((10, 5, 8))
    min_grididx = (1, 1, 1)
    quantics = QuanticsGrids.grididx_to_quantics(grid, min_grididx)
    @test all(q -> q == 1, quantics)
    @test QuanticsGrids.quantics_to_grididx(grid, quantics) == min_grididx

    # Test with maximum valid grididx
    max_grididx = (2^10, 2^5, 2^8)
    quantics = QuanticsGrids.grididx_to_quantics(grid, max_grididx)
    @test all(q -> q == 2, quantics)
    @test QuanticsGrids.quantics_to_grididx(grid, quantics) == max_grididx
end

@testitem "challenging tests - mixed bases" begin
    # Test with base 3
    grid = QuanticsGrids.NewDiscretizedGrid((4, 6, 3); base=3)
    for _ in 1:50
        quantics = rand(1:3, length(grid))
        grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
        @test QuanticsGrids.grididx_to_quantics(grid, grididx) == quantics
        @test all(1 .<= grididx .<= (3 .^ grid.Rs))
    end

    # Test with base 5
    grid = QuanticsGrids.NewDiscretizedGrid((3, 4); base=5)
    for _ in 1:50
        grididx = (rand(1:5^3), rand(1:5^4))
        quantics = QuanticsGrids.grididx_to_quantics(grid, grididx)
        @test QuanticsGrids.quantics_to_grididx(grid, quantics) == grididx
        @test all(1 .<= quantics .<= 5)
    end
end

@testitem "challenging tests - complex fused indices" begin
    # Multiple variables fused in single sites
    grid = QuanticsGrids.NewDiscretizedGrid(
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
    @test QuanticsGrids.quantics_to_grididx(grid, [8, 1, 4, 2, 1]) == (6, 4, 7)
    @test QuanticsGrids.grididx_to_quantics(grid, (6, 4, 7)) == [8, 1, 4, 2, 1]

    # Test random values
    for _ in 1:100
        quantics = [
            rand(1:8),   # site 1: base^3 = 8 possibilities
            rand(1:2),   # site 2: base = 2 possibilities  
            rand(1:4),   # site 3: base^2 = 4 possibilities
            rand(1:2),   # site 4: base = 2 possibilities
            rand(1:2)    # site 5: base = 2 possibilities
        ]
        grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
        @test QuanticsGrids.grididx_to_quantics(grid, grididx) == quantics
    end
end

@testitem "challenging tests - asymmetric grids" begin
    # Very asymmetric grid dimensions
    grid = QuanticsGrids.NewDiscretizedGrid((20, 0, 15, 3))

    for _ in 1:50
        # Test boundary conditions
        quantics = rand(1:2, length(grid))
        grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
        @test QuanticsGrids.grididx_to_quantics(grid, grididx) == quantics

        # Verify grid index bounds
        @test 1 <= grididx[1] <= 2^20
        @test grididx[2] == 1  # R=0 means only one possible value
        @test 1 <= grididx[3] <= 2^15
        @test 1 <= grididx[4] <= 2^3
    end
end

@testitem "challenging tests - single dimension edge cases" begin
    # 1D grid with large R
    grid = QuanticsGrids.NewDiscretizedGrid((25,))

    # Test extremes
    min_quantics = ones(Int, 25)
    max_quantics = fill(2, 25)

    @test QuanticsGrids.quantics_to_grididx(grid, min_quantics) == (1,)
    @test QuanticsGrids.quantics_to_grididx(grid, max_quantics) == (2^25,)
    @test QuanticsGrids.grididx_to_quantics(grid, (1,)) == min_quantics
    @test QuanticsGrids.grididx_to_quantics(grid, (2^25,)) == max_quantics

    # Test middle values
    for _ in 1:20
        quantics = rand(1:2, 25)
        grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
        @test QuanticsGrids.grididx_to_quantics(grid, grididx) == quantics
    end
end

@testitem "challenging tests - high dimensional grids" begin
    # 8D grid with moderate R values
    grid = QuanticsGrids.NewDiscretizedGrid(ntuple(i -> 4 + (i % 3), 8); base=3)

    for _ in 1:30
        quantics = rand(1:3, length(grid))
        grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
        @test QuanticsGrids.grididx_to_quantics(grid, grididx) == quantics

        # Verify all grid indices are within bounds
        for d in 1:8
            @test 1 <= grididx[d] <= 3^grid.Rs[d]
        end
    end
end

@testitem "challenging tests - stress test with complex patterns" begin
    # Grid with intentionally complex index table structure
    grid = QuanticsGrids.NewDiscretizedGrid(
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

        grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
        recovered = QuanticsGrids.grididx_to_quantics(grid, grididx)
        @test recovered == quantics

        # Verify grid indices are reasonable
        @test all(1 .<= grididx .<= (3 .^ grid.Rs))
    end
end

# @testset "comprehensive quantics_to_grididx tests" begin
#     @testset "simple 1D cases" begin
#         grid = QuanticsGrids.NewDiscretizedGrid((4,))
#         @test QuanticsGrids.quantics_to_grididx(grid, [1, 1, 1, 1]) == (1,)
#         @test QuanticsGrids.quantics_to_grididx(grid, [2, 1, 1, 1]) == (9,)
#         @test QuanticsGrids.quantics_to_grididx(grid, [1, 2, 1, 1]) == (5,)
#         @test QuanticsGrids.quantics_to_grididx(grid, [2, 2, 2, 2]) == (16,)
#     end

#     @testset "2D rectangular grids" begin
#         grid = QuanticsGrids.NewDiscretizedGrid((2, 3))
#         # Test all possible combinations
#         @test QuanticsGrids.quantics_to_grididx(grid, [1, 1, 1, 1, 1]) == (1, 1)
#         @test QuanticsGrids.quantics_to_grididx(grid, [2, 2, 1, 1, 1]) == (3, 5)
#         @test QuanticsGrids.quantics_to_grididx(grid, [1, 1, 2, 1, 1]) == (2, 1)
#         @test QuanticsGrids.quantics_to_grididx(grid, [2, 2, 2, 2, 2]) == (4, 8)
#     end

#     @testset "different bases" begin
#         # Base 3 tests
#         grid = QuanticsGrids.NewDiscretizedGrid((2, 2); base=3)
#         @test QuanticsGrids.quantics_to_grididx(grid, [1, 1, 1, 1]) == (1, 1)
#         @test QuanticsGrids.quantics_to_grididx(grid, [3, 1, 1, 1]) == (7, 1)
#         @test QuanticsGrids.quantics_to_grididx(grid, [1, 3, 1, 1]) == (1, 7)
#         @test QuanticsGrids.quantics_to_grididx(grid, [3, 3, 3, 3]) == (9, 9)

#         # Base 4 tests
#         grid = QuanticsGrids.NewDiscretizedGrid((1, 2); base=4)
#         @test QuanticsGrids.quantics_to_grididx(grid, [1, 1, 1]) == (1, 1)
#         @test QuanticsGrids.quantics_to_grididx(grid, [4, 4, 4]) == (4, 16)
#     end

#     @testset "fused indices complexity" begin
#         # Test with multiple fused combinations
#         grid = QuanticsGrids.NewDiscretizedGrid(
#             (:w, :x, :y, :z),
#             [
#                 [(:w, 2), (:x, 1)],     # site 1: 2 vars
#                 [(:y, 3)],              # site 2: 1 var
#                 [(:w, 1), (:z, 2)],     # site 3: 2 vars
#                 [(:x, 2), (:y, 2), (:z, 1)], # site 4: 3 vars
#                 [(:y, 1)]               # site 5: 1 var
#             ]
#         )

#         # Test specific patterns
#         @test QuanticsGrids.quantics_to_grididx(grid, [1, 1, 1, 1, 1]) == (1, 1, 1, 1)
#         @test QuanticsGrids.quantics_to_grididx(grid, [4, 2, 4, 8, 2]) == (4, 4, 8, 4)

#         # Test boundary values
#         quantics_min = [1, 1, 1, 1, 1]
#         quantics_max = [4, 2, 4, 8, 2]
#         grididx_min = QuanticsGrids.quantics_to_grididx(grid, quantics_min)
#         grididx_max = QuanticsGrids.quantics_to_grididx(grid, quantics_max)
#         @test all(grididx_min .== (1, 1, 1, 1))
#         @test all(grididx_max .== 2 .^ grid.Rs)
#     end

#     @testset "pathological cases" begin
#         # Very unbalanced Rs
#         grid = QuanticsGrids.NewDiscretizedGrid((1, 10, 1, 5))
#         quantics = [1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1]
#         grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
#         @test grididx[1] == 1  # R=1 means only one value
#         @test grididx[3] == 1  # R=1 means only one value
#         @test 1 <= grididx[2] <= 2^10
#         @test 1 <= grididx[4] <= 2^5
#     end
# end

# @testset "comprehensive grididx_to_quantics tests" begin
#     @testset "simple 1D reconstruction" begin
#         grid = QuanticsGrids.NewDiscretizedGrid((3,))
#         @test QuanticsGrids.grididx_to_quantics(grid, (1,)) == [1, 1, 1]
#         @test QuanticsGrids.grididx_to_quantics(grid, (8,)) == [2, 2, 2]
#         @test QuanticsGrids.grididx_to_quantics(grid, (5,)) == [2, 1, 1]
#         @test QuanticsGrids.grididx_to_quantics(grid, (3,)) == [1, 2, 1]
#     end

#     @testset "2D reconstruction" begin
#         grid = QuanticsGrids.NewDiscretizedGrid((2, 3))
#         @test QuanticsGrids.grididx_to_quantics(grid, (1, 1)) == [1, 1, 1, 1, 1]
#         @test QuanticsGrids.grididx_to_quantics(grid, (4, 8)) == [2, 2, 2, 2, 2]
#         @test QuanticsGrids.grididx_to_quantics(grid, (2, 3)) == [1, 1, 2, 2, 1]
#     end

#     @testset "base conversion accuracy" begin
#         # Test with base 5
#         grid = QuanticsGrids.NewDiscretizedGrid((3, 2); base=5)
#         test_cases = [
#             ((1, 1), [1, 1, 1, 1, 1]),
#             ((125, 25), [5, 5, 5, 5, 5]),
#             ((63, 13), [3, 3, 3, 3, 3]),
#             ((31, 6), [2, 2, 2, 1, 1])
#         ]

#         for (grididx, expected_quantics) in test_cases
#             @test QuanticsGrids.grididx_to_quantics(grid, grididx) == expected_quantics
#         end
#     end

#     @testset "complex fused reconstruction" begin
#         # Highly complex fused structure
#         grid = QuanticsGrids.NewDiscretizedGrid(
#             (:p, :q, :r, :s, :t),
#             [
#                 [(:p, 3)],                      # site 1
#                 [(:q, 2), (:r, 4)],            # site 2: 2 vars
#                 [(:s, 1)],                      # site 3
#                 [(:p, 2), (:t, 3)],            # site 4: 2 vars
#                 [(:q, 1), (:r, 3), (:s, 2)],  # site 5: 3 vars
#                 [(:p, 1)],                      # site 6
#                 [(:r, 2)],                      # site 7
#                 [(:t, 2)],                      # site 8
#                 [(:r, 1), (:s, 3), (:t, 1)]   # site 9: 3 vars
#             ];
#             base=3
#         )

#         # Test specific known reconstructions
#         @test QuanticsGrids.grididx_to_quantics(grid, (1, 1, 1, 1, 1)) == [1, 1, 1, 1, 1, 1, 1, 1, 1]
#         @test QuanticsGrids.grididx_to_quantics(grid, (27, 9, 27, 27, 27)) == [3, 9, 3, 9, 27, 3, 3, 3, 9]

#         # Test some intermediate values
#         grididx = (14, 5, 3, 7, 15)
#         quantics = QuanticsGrids.grididx_to_quantics(grid, grididx)
#         @test length(quantics) == 9
#         @test all(1 .<= quantics .<= QuanticsGrids.sitedim.(Ref(grid), 1:length(grid)))
#     end

#     @testset "boundary and edge value reconstruction" begin
#         grid = QuanticsGrids.NewDiscretizedGrid((8, 4, 6); base=4)

#         # Test minimum values
#         min_grididx = (1, 1, 1)
#         min_quantics = QuanticsGrids.grididx_to_quantics(grid, min_grididx)
#         @test all(q -> q == 1, min_quantics)

#         # Test maximum values
#         max_grididx = (4^8, 4^4, 4^6)
#         max_quantics = QuanticsGrids.grididx_to_quantics(grid, max_grididx)
#         @test all(q -> q == 4, max_quantics)

#         # Test power-of-base values
#         for d in 1:3
#             for exp in 1:grid.Rs[d]-1
#                 grididx = ntuple(i -> i == d ? 4^exp + 1 : 1, 3)
#                 quantics = QuanticsGrids.grididx_to_quantics(grid, grididx)
#                 @test length(quantics) == sum(grid.Rs)
#                 @test all(1 .<= quantics .<= 4)
#             end
#         end
#     end
# end

# @testset "round-trip consistency under stress" begin
#     @testset "massive random testing" begin
#         # Large scale random testing with different configurations
#         configs = [
#             ((5, 7, 3), 2),
#             ((4, 4, 4, 4), 3),
#             ((10, 2, 8), 2),
#             ((3, 6, 4, 5), 4),
#             ((12, 1, 15), 2)
#         ]

#         for (Rs, base) in configs
#             grid = QuanticsGrids.NewDiscretizedGrid(Rs; base)

#             # Test quantics -> grididx -> quantics
#             for _ in 1:100
#                 quantics = rand(1:base, length(grid))
#                 grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
#                 recovered = QuanticsGrids.grididx_to_quantics(grid, grididx)
#                 @test recovered == quantics
#             end

#             # Test grididx -> quantics -> grididx
#             for _ in 1:100
#                 grididx = ntuple(d -> rand(1:base^Rs[d]), length(Rs))
#                 quantics = QuanticsGrids.grididx_to_quantics(grid, grididx)
#                 recovered = QuanticsGrids.quantics_to_grididx(grid, quantics)
#                 @test recovered == grididx
#             end
#         end
#     end

#     @testset "systematic boundary testing" begin
#         grid = QuanticsGrids.NewDiscretizedGrid((4, 5, 3); base=3)

#         # Test all corner cases systematically
#         for d1 in [1, 3^4]
#             for d2 in [1, 3^5]
#                 for d3 in [1, 3^3]
#                     grididx = (d1, d2, d3)
#                     quantics = QuanticsGrids.grididx_to_quantics(grid, grididx)
#                     recovered = QuanticsGrids.quantics_to_grididx(grid, quantics)
#                     @test recovered == grididx
#                 end
#             end
#         end
#     end

#     @testset "fused indices stress test" begin
#         # Extremely complex fused structure
#         grid = QuanticsGrids.NewDiscretizedGrid(
#             (:u, :v, :w, :x, :y, :z),
#             [
#                 [(:u, 4), (:v, 3), (:w, 2)],           # site 1: 3 vars
#                 [(:x, 5)],                               # site 2: 1 var
#                 [(:y, 4), (:z, 4)],                     # site 3: 2 vars
#                 [(:u, 3), (:v, 2)],                     # site 4: 2 vars
#                 [(:w, 1), (:x, 4), (:y, 3), (:z, 3)],  # site 5: 4 vars
#                 [(:u, 2)],                               # site 6: 1 var
#                 [(:v, 1), (:x, 3)],                     # site 7: 2 vars
#                 [(:w, 3), (:y, 2)],                     # site 8: 2 vars
#                 [(:x, 2)],                               # site 9: 1 var
#                 [(:u, 1), (:z, 2)],                     # site 10: 2 vars
#                 [(:y, 1)],                               # site 11: 1 var
#                 [(:z, 1)],                               # site 12: 1 var
#                 [(:x, 1), (:w, 4)]                      # site 13: 2 vars
#             ];
#             base=5
#         )

#         # Intensive round-trip testing
#         for _ in 1:500
#             # Generate valid quantics
#             quantics = [
#                 rand(1:125),  # site 1: 5^3
#                 rand(1:5),    # site 2: 5^1
#                 rand(1:25),   # site 3: 5^2
#                 rand(1:25),   # site 4: 5^2
#                 rand(1:625),  # site 5: 5^4
#                 rand(1:5),    # site 6: 5^1
#                 rand(1:25),   # site 7: 5^2
#                 rand(1:25),   # site 8: 5^2
#                 rand(1:5),    # site 9: 5^1
#                 rand(1:25),   # site 10: 5^2
#                 rand(1:5),    # site 11: 5^1
#                 rand(1:5),    # site 12: 5^1
#                 rand(1:25)    # site 13: 5^2
#             ]

#             grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
#             recovered = QuanticsGrids.grididx_to_quantics(grid, grididx)
#             @test recovered == quantics

#             # Verify grididx bounds
#             @test all(1 .<= grididx .<= (5 .^ grid.Rs))
#         end
#     end
# end

# @testset "performance and correctness under extreme conditions" begin
#     @testset "large dimension count" begin
#         # Test with many dimensions but small R values
#         Rs = ntuple(i -> 2 + (i % 3), 12)  # 12 dimensions
#         grid = QuanticsGrids.NewDiscretizedGrid(Rs; base=2)

#         for _ in 1:50
#             quantics = rand(1:2, length(grid))
#             grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
#             recovered = QuanticsGrids.grididx_to_quantics(grid, grididx)
#             @test recovered == quantics
#             @test length(grididx) == 12
#         end
#     end

#     @testset "large R values with small dimension count" begin
#         # Test with few dimensions but large R values
#         grid = QuanticsGrids.NewDiscretizedGrid((18, 16); base=2)

#         # Test specific patterns
#         patterns = [
#             ones(Int, length(grid)),
#             fill(2, length(grid)),
#             [rand(1:2) for _ in 1:length(grid)]
#         ]

#         for pattern in patterns
#             grididx = QuanticsGrids.quantics_to_grididx(grid, pattern)
#             recovered = QuanticsGrids.grididx_to_quantics(grid, grididx)
#             @test recovered == pattern
#         end
#     end

#     @testset "mixed high-low R combinations" begin
#         grid = QuanticsGrids.NewDiscretizedGrid((0, 15, 0, 12, 0, 8); base=3)

#         for _ in 1:100
#             quantics = rand(1:3, length(grid))
#             grididx = QuanticsGrids.quantics_to_grididx(grid, quantics)
#             recovered = QuanticsGrids.grididx_to_quantics(grid, grididx)
#             @test recovered == quantics

#             # Check that dimensions with R=0 always give grididx=1
#             @test grididx[1] == 1
#             @test grididx[3] == 1
#             @test grididx[5] == 1
#         end
#     end
# end
