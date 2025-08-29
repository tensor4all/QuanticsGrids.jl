@testitem "quantics.jl" begin
    import QuanticsGrids: fuse_dimensions, unfuse_dimensions
    import QuanticsGrids: quantics_to_index_fused, index_to_quantics
    import QuanticsGrids: interleave_dimensions, deinterleave_dimensions
    import QuanticsGrids: serialize_dimensions, deserialize_dimensions
    import QuanticsGrids: meander_dimensions, demeander_dimensions

    @testset "quantics representation" begin
        @testset "fuse_dimensions" begin
            x = [1, 1]
            y = [2, 2]
            fused = fuse_dimensions(x, y)
            @test unfuse_dimensions(fused, 2)[1] == x
            @test unfuse_dimensions(fused, 2)[2] == y
        end

        @testset "quantics to index, 1d" begin
            B = 2
            @test quantics_to_index_fused([1, 1, 1, 1]; base = B) == (1,)
            @test quantics_to_index_fused([1, 1, 1, 2]; base = B) == (2,)
            @test quantics_to_index_fused([1, 1, 2, 1]; base = B) == (3,)
            @test quantics_to_index_fused([1, 1, 2, 2]; base = B) == (4,)
            @test quantics_to_index_fused([2, 1, 1, 1]; base = B) == (9,)
        end

        @testset "quantics to index, 1d (general power base)" for B in [2, 4]
            d = 1
            R = 4 # number of digits

            index_reconst = Int[]
            for index = 1:B^R
                digitlist_ = QuanticsGrids.index_to_quantics(index; numdigits = R, base = B)
                push!(
                    index_reconst,
                    only(quantics_to_index_fused(digitlist_; base = B, dims = Val(d))),
                )
            end

            @test collect(1:B^R) == index_reconst
        end

        @testset "quantics to index, 2d, base=3" begin
            base = 3
            dim = 2

            # X_i = Fused quantics index at i (1 <= i <= R)
            # x_i = quantics index for the first variable at i (1 <= i <= R)
            # y_i = quantics index for the second variable at i (1 <= i <= R)
            #
            # X_i = (x_i-1) + (base) * (y_i-1) + 1 (column major)
            @test quantics_to_index_fused([1, 1]; base = base, dims = Val(dim)) == (1, 1)
            @test quantics_to_index_fused([1, 2]; base = base, dims = Val(dim)) == (2, 1)
            @test quantics_to_index_fused([1, 3]; base = base, dims = Val(dim)) == (3, 1)
            @test quantics_to_index_fused([1, 4]; base = base, dims = Val(dim)) == (1, 2)
            @test quantics_to_index_fused([1, 5]; base = base, dims = Val(dim)) == (2, 2)
            @test quantics_to_index_fused([1, 6]; base = base, dims = Val(dim)) == (3, 2)
            @test quantics_to_index_fused([1, 7]; base = base, dims = Val(dim)) == (1, 3)
            @test quantics_to_index_fused([1, 8]; base = base, dims = Val(dim)) == (2, 3)
            @test quantics_to_index_fused([1, 9]; base = base, dims = Val(dim)) == (3, 3)
            @test quantics_to_index_fused([2, 1]; base = base, dims = Val(dim)) == (4, 1)
        end

        @testset "quantics back-and-forth, 2d" begin
            base = 2
            dim = 2
            R = 2

            for j = 1:base^R, i = 1:base^R
                index = (i, j)

                digitlist1 = Vector{Int}(undef, R)
                QuanticsGrids.index_to_quantics_fused!(digitlist1, index; base = base)
                digitlist2 =
                    QuanticsGrids.index_to_quantics_fused(index; numdigits = R, base = base)
                @test digitlist1 == digitlist2

                index_reconst = QuanticsGrids.quantics_to_index_fused(
                    digitlist1;
                    base = base,
                    dims = Val(dim),
                )
                @test index == index_reconst
            end
        end

        @testset "interleave dimensions" begin
            @test [1, 1, 1, 1] == interleave_dimensions([1, 1, 1, 1])
            @test [1, 1, 1, 1, 1, 1, 1, 1] ==
                  interleave_dimensions([1, 1, 1, 1], [1, 1, 1, 1])
            @test [1, 2, 1, 3, 1, 4, 1, 5] ==
                  interleave_dimensions([1, 1, 1, 1], [2, 3, 4, 5])
            @test [1, 2, 11, 1, 3, 12, 1, 4, 13, 1, 5, 14] ==
                  interleave_dimensions([1, 1, 1, 1], [2, 3, 4, 5], [11, 12, 13, 14])
        end

        @testset "deinterleave dimensions" begin
            @test deinterleave_dimensions([1, 1, 1, 1], 1) == [[1, 1, 1, 1]]
            @test deinterleave_dimensions([1, 1, 1, 1], 2) == [[1, 1], [1, 1]]
            @test deinterleave_dimensions([1, 2, 1, 3, 1, 4, 1, 5], 2) ==
                  [[1, 1, 1, 1], [2, 3, 4, 5]]
            @test deinterleave_dimensions([1, 2, 11, 1, 3, 12, 1, 4, 13, 1, 5, 14], 3) ==
                  [[1, 1, 1, 1], [2, 3, 4, 5], [11, 12, 13, 14]]
        end

        @testset "serial dimensions" begin
            @test [1, 1, 1, 1] == serialize_dimensions([1, 1, 1, 1])
            @test [1, 1, 1, 1, 1, 1, 1, 1] ==
                  serialize_dimensions([1, 1, 1, 1], [1, 1, 1, 1])
            @test [1, 1, 1, 1, 2, 3, 4, 5] ==
                  serialize_dimensions([1, 1, 1, 1], [2, 3, 4, 5])
            @test [1, 1, 1, 1, 2, 3, 4, 5, 11, 12, 13, 14] ==
                  serialize_dimensions([1, 1, 1, 1], [2, 3, 4, 5], [11, 12, 13, 14])
        end

        @testset "deserial dimensions" begin
            @test deserialize_dimensions([1, 1, 1, 1], 1) == [[1, 1, 1, 1]]
            @test deserialize_dimensions([1, 1, 1, 1], 2) == [[1, 1], [1, 1]]
            @test deserialize_dimensions([1, 2, 1, 3, 1, 4, 1, 5], 2) ==
                  [[1, 2, 1, 3], [1, 4, 1, 5]]
            @test deserialize_dimensions([1, 2, 11, 1, 3, 12, 1, 4, 13, 1, 5, 14], 3) ==
                  [[1, 2, 11, 1], [3, 12, 1, 4], [13, 1, 5, 14]]
        end

        @testset "meander dimensions" begin
            @test [1, 1, 1, 1] == meander_dimensions([1, 1, 1, 1])
            @test [1, 1, 1, 1, 1, 1, 1, 1] ==
                  meander_dimensions([1, 1, 1, 1], [1, 1, 1, 1])
            @test [1, 1, 1, 1, 2, 3, 4, 5] ==
                  meander_dimensions([1, 1, 1, 1], [2, 3, 4, 5])
            @test [1, 1, 1, 1, 2, 3, 4, 5, 11, 12, 13, 14] ==
                  meander_dimensions([1, 1, 1, 1], [2, 3, 4, 5], [14, 13, 12, 11])
        end

        @testset "demeander dimensions" begin
            @test demeander_dimensions([1, 1, 1, 1], 1) == [[1, 1, 1, 1]]
            @test demeander_dimensions([1, 1, 1, 1], 2) == [[1, 1], [1, 1]]
            @test demeander_dimensions([1, 2, 1, 3, 1, 4, 1, 5], 2) ==
                  [[3, 1, 2, 1], [1, 4, 1, 5]]
            @test demeander_dimensions([1, 2, 11, 1, 3, 12, 1, 4, 13, 1, 5, 14], 3) ==
                  [[1, 11, 2, 1], [3, 12, 1, 4], [14, 5, 1, 13]]
        end
    end
end
