@testitem "Code quality (Aqua.jl)" begin
    using Aqua

    Aqua.test_all(QuanticsGrids; deps_compat=false)
end

@testitem "Code linting (JET.jl)" begin
    using JET

    if VERSION >= v"1.10"
        JET.test_package(QuanticsGrids; target_modules=(QuanticsGrids,))
    end
end
