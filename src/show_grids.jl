
# ============================================================================
# Display/show methods
# ============================================================================

function Base.show(io::IO, ::MIME"text/plain", g::DiscretizedGrid{D}) where D
    print(io, "DiscretizedGrid{$D}")

    # Grid resolution and total points
    total_points = prod(grid_base(g) .^ grid_Rs(g))
    if D == 1
        print(io, " with $(grid_base(g)^grid_Rs(g)[1]) grid points")
    else
        print(io, " with $(join(grid_base(g) .^ grid_Rs(g), "×")) = $total_points grid points")
    end

    # Variable names (if meaningful)
    if any(name -> !startswith(string(name), r"^\d+$"), grid_variablenames(g))
        var_str = join(grid_variablenames(g), ", ")
        print(io, "\n├─ Variables: ($var_str)")
    end

    # Resolution per dimension
    if D == 1
        print(io, "\n├─ Resolution: $(grid_Rs(g)[1]) bits")
    else
        res_str = join(["$(grid_variablenames(g)[i]): $(grid_Rs(g)[i])" for i in 1:D], ", ")
        print(io, "\n├─ Resolutions: ($res_str)")
    end

    # Bounds (only show if not default unit interval/square/cube)
    default_lower = default_lower_bound(Val(D))
    default_upper = default_upper_bound(Val(D))
    if lower_bound(g) != default_lower || any(abs.(upper_bound(g) .- default_upper) .> 1e-10)
        if D == 1
            print(io, "\n├─ Domain: [$(lower_bound(g)[1]), $(upper_bound(g)[1]))")
        else
            bounds_str = join(["[$(lower_bound(g)[i]), $(upper_bound(g)[i]))" for i in 1:D], " × ")
            print(io, "\n├─ Domain: $bounds_str")
        end

        # Grid spacing
        step_vals = grid_step(g)
        if D == 1
            print(io, "\n├─ Grid spacing: $(step_vals)")
        else
            step_str = join(["Δ$(grid_variablenames(g)[i]) = $(step_vals[i])" for i in 1:D], ", ")
            print(io, "\n├─ Grid spacing: ($step_str)")
        end
    else
        # For unit domain, show appropriate unit description
        unit_domain_str = if D == 1
            "unit interval [0, 1)"
        elseif D == 2
            "unit square [0, 1)²"
        elseif D == 3
            "unit cube [0, 1)³"
        else
            "unit hypercube [0, 1)^$D"
        end
        print(io, "\n├─ Domain: $unit_domain_str")

        # Grid spacing
        step_vals = grid_step(g)
        if D == 1
            print(io, "\n├─ Grid spacing: $(step_vals)")
        else
            step_str = join(["Δ$(grid_variablenames(g)[i]) = $(step_vals[i])" for i in 1:D], ", ")
            print(io, "\n├─ Grid spacing: ($step_str)")
        end
    end

    # Base (only show if not binary)
    if grid_base(g) != 2
        print(io, "\n├─ Base: $(grid_base(g))")
    end

    # Tensor structure summary
    num_sites = length(grid_indextable(g))
    site_dims = [grid_base(g)^length(site) for site in grid_indextable(g)]
    max_bond_dim = maximum(site_dims)

    print(io, "\n└─ Tensor train: $num_sites sites")
    if all(d -> d == site_dims[1], site_dims)
        print(io, " (uniform dimension $(site_dims[1]))")
    else
        print(io, " (dimensions: $(join(site_dims, "-")))")
    end
end

function Base.show(io::IO, g::DiscretizedGrid{D}) where D
    total_points = prod(grid_base(g) .^ grid_Rs(g))
    if D == 1
        print(io, "DiscretizedGrid{$D}($(grid_base(g)^grid_Rs(g)[1]) points)")
    else
        print(io, "DiscretizedGrid{$D}($(join(grid_base(g) .^ grid_Rs(g), "×")) points)")
    end
end
