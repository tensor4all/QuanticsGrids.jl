@testitem "finite_difference.jl" begin
    using Test
    import QuanticsGrids
    using QuanticsGrids: central_difference_periodic

    # Grid points as a Vector (the package works with ranges, but a Vector keeps
    # the arithmetic below free of range optimisations).
    gridpoints(g) = collect(QuanticsGrids.grid_origcoords(g, 1))

    @testset "type and size" begin
        g = QuanticsGrids.DiscretizedGrid{1}(4, -0.5, 0.5)
        for order in (1, 2)
            D = central_difference_periodic(g; order)
            @test D isa Matrix{Float64}
            @test size(D) == (16, 16)
            # the derivative annihilates constants
            @test all(iszero, D * ones(16))
        end
        @test central_difference_periodic(g) == central_difference_periodic(g; order=1)
    end

    @testset "overlapping stencils on tiny grids" begin
        # n = 1 and n = 2 are the cases where both neighbours are the same
        # grid point, so the contributions have to add up instead of overwrite.
        g1 = QuanticsGrids.DiscretizedGrid{1}(0, -1.0, 1.0)
        @test central_difference_periodic(g1) == zeros(1, 1)
        @test central_difference_periodic(g1; order=2) == zeros(1, 1)
        g2 = QuanticsGrids.DiscretizedGrid{1}(1, -1.0, 1.0)
        Δ2 = QuanticsGrids.grid_step(g2)
        @test central_difference_periodic(g2) == zeros(2, 2)
        @test central_difference_periodic(g2; order=2) == [-2 2; 2 -2] ./ Δ2^2
    end

    @testset "stencil coefficients come from grid_step (d = 1 and d = 4)" begin
        for (lower, upper) in ((-0.5, 0.5), (-2.0, 2.0)), base in (2, 3)
            g = QuanticsGrids.DiscretizedGrid{1}(3, lower, upper; base)
            n = base^3
            Δ = QuanticsGrids.grid_step(g)
            @test Δ ≈ (upper - lower) / n
            @test n == length(gridpoints(g))
            D1 = central_difference_periodic(g; order=1)
            D2 = central_difference_periodic(g; order=2)
            # An implementation that hard-codes Δ = base^-R instead of using
            # grid_step(g) would give coefficients that are too large by d.
            @test D1[2, 3] == 1 / (2Δ)
            @test D1[2, 1] == -1 / (2Δ)
            @test D2[2, 2] == -2 / Δ^2
            @test D2[2, 3] == 1 / Δ^2
            @test D2[2, 1] == 1 / Δ^2
        end
    end

    @testset "linear function a*v + b (d = 1 and d = 4)" begin
        for (lower, upper) in ((-0.5, 0.5), (-2.0, 2.0))
            g = QuanticsGrids.DiscretizedGrid{1}(4, lower, upper)
            n = 16
            Δ = QuanticsGrids.grid_step(g)
            d = upper - lower
            x = gridpoints(g)
            a, b = 3.0, -1.0
            f = a .* x .+ b
            D1 = central_difference_periodic(g; order=1)
            D2 = central_difference_periodic(g; order=2)

            # A linear function is not periodic, so the two rows that wrap
            # around see the jump across the boundary and are excluded here;
            # they are checked below. These fixtures are exact in Float64, so
            # the derivative is required to be exactly `a` and exactly 0.
            @test all(==(a), (D1 * f)[2:(n-1)])
            @test all(iszero, (D2 * f)[2:(n-1)])

            # The wrap rows differ from `a` by a *periodicity* effect: the left
            # neighbour of x_1 is x_n = x_1 - Δ + d, not x_1 - Δ. Both wrapped
            # rows see the same jump 2Δ - d of the linear function.
            @test (D1 * f)[1] ≈ a * (2Δ - d) / (2Δ)
            @test (D1 * f)[n] ≈ a * (2Δ - d) / (2Δ)
            @test (D2 * f)[1] ≈ a * d / Δ^2
            @test (D2 * f)[n] ≈ -a * d / Δ^2
        end
    end

    @testset "periodic functions: second-order convergence" begin
        # Periodic sin/cos on grids of different physical size; the analytic
        # derivatives use the physical wavenumber 2π/d, so a wrong Δ shows up
        # both in the prefactor and in the convergence.
        for (lower, upper) in ((-0.5, 0.5), (-2.0, 2.0))
            d = upper - lower
            k = 2π / d
            errs1 = Float64[]
            errs2 = Float64[]
            for R in 4:8
                g = QuanticsGrids.DiscretizedGrid{1}(R, lower, upper)
                x = gridpoints(g)
                for (f, df, d2f) in (
                    (sin.(k .* x), k .* cos.(k .* x), -k^2 .* sin.(k .* x)),
                    (cos.(k .* x), -k .* sin.(k .* x), -k^2 .* cos.(k .* x)),
                )
                    push!(errs1, maximum(abs, central_difference_periodic(g; order=1) * f .- df))
                    push!(errs2, maximum(abs, central_difference_periodic(g; order=2) * f .- d2f))
                end
            end
            # error = O(Δ^2): halving Δ quarters the error (same function, one
            # more resolution; sin and cos occupy consecutive entries)
            for errs in (errs1, errs2)
                @test all(>(0), errs)
                for j in 1:2:(length(errs)-2)
                    @test errs[j] / errs[j+2] ≈ 4 rtol=0.05
                end
            end
        end
    end

    @testset "periodic wrap (i = 1 and i = n)" begin
        g = QuanticsGrids.DiscretizedGrid{1}(3, -2.0, 2.0)
        n = 8
        Δ = QuanticsGrids.grid_step(g)
        D1 = central_difference_periodic(g; order=1)
        D2 = central_difference_periodic(g; order=2)

        # The operator is circulant (row 1 is the reference row).
        for i in 2:n, j in 1:n
            @test D1[i, j] == D1[1, mod1(j - i + 1, n)]
            @test D2[i, j] == D2[1, mod1(j - i + 1, n)]
        end

        # A delta at the last grid point is felt at i = 1 and i = n - 1.
        e = zeros(n)
        e[n] = 1.0
        v = zeros(n)
        v[1] = -1 / (2Δ)
        v[n-1] = 1 / (2Δ)
        @test D1 * e ≈ v
        w = zeros(n)
        w[1] = 1 / Δ^2
        w[n-1] = 1 / Δ^2
        w[n] = -2 / Δ^2
        @test D2 * e ≈ w

        # ... and one at the first grid point at i = 2 and i = n.
        e = zeros(n)
        e[1] = 1.0
        v = zeros(n)
        v[2] = -1 / (2Δ)
        v[n] = 1 / (2Δ)
        @test D1 * e ≈ v
        w = zeros(n)
        w[1] = -2 / Δ^2
        w[2] = 1 / Δ^2
        w[n] = 1 / Δ^2
        @test D2 * e ≈ w
    end

    @testset "argument validation" begin
        g = QuanticsGrids.DiscretizedGrid{1}(3, -1.0, 1.0)
        @test_throws ArgumentError central_difference_periodic(g; order=3)
        @test_throws ArgumentError central_difference_periodic(g; order=0)
        g_endpoint = QuanticsGrids.DiscretizedGrid{1}(3, -1.0, 1.0; includeendpoint=true)
        @test_throws ArgumentError central_difference_periodic(g_endpoint)
    end
end
