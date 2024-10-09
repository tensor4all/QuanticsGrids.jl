@testitem "Code quality (Aqua.jl)" begin
    using Aqua

    import QuanticsGrids

    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(QuanticsGrids; unbound_args = false, deps_compat = false)
    end

end

@testitem "Code linting (JET.jl)" begin
    using Test
    using JET

    import QuanticsGrids

    if VERSION >= v"1.9"
        @testset "Code linting (JET.jl)" begin
            JET.test_package(QuanticsGrids; target_defined_modules = true)
        end
    end
end
