export NewDiscretizedGrid

"""
    NewDiscretizedGrid{D}

A discretized grid structure for D-dimensional grids with variable resolution,
supporting efficient conversion between quantics, grid indices, and original coordinates.
A `NewDiscretizedGrid` instance is intended to undergird a quantics tensor train
with a specific index structure, as defined in the `indextable` field.
For example, say indextable is `[[(:a, 1), (:b, 2)], [(:a, 2)], [(:b, 1), (:a, 3)]]`,
then the corresponding tensor train has 3 tensor cores:

      a_1 b_2          a_2          b_1 a_3
       |   |            |            |   |
    ┌──┴───┴──┐    ┌────┴────┐    ┌──┴───┴──┐
    │         │    │         │    │         │
    │         │────│         │────│         │
    │         │    │         │    │         │
    └─────────┘    └─────────┘    └─────────┘

This object may be constructed with
```julia-repl
julia> grid = NewDiscretizedGrid((:a, :b), [[(:a, 1), (:b, 2)], [(:a, 2)], [(:b, 1), (:a, 3)]])
NewDiscretizedGrid{2} with 8×4 = 32 grid points
├─ Variables: (a, b)
├─ Resolutions: (a: 3, b: 2)
├─ Domain: unit square [0, 1)²
├─ Grid spacing: (Δa = 0.125, Δb = 0.25)
└─ Tensor train: 3 sites (dimensions: 4-2-4)
```
and represents a 2^3 x 2^2 discretization of the unit square in the 2D plane (the x mark grid points):

       1.0  ┌───────────────────────────────┐
            │                               │
            │                               │
       0.75 x   x   x   x   x   x   x   x   │
            │                               │
            │                               │
    b  0.5  x   x   x   x   x   x   x   x   │
            │                               │
            │                               │
       0.25 x   x   x   x   x   x   x   x   │
            │                               │
            │                               │
       0.0  x───x───x───x───x───x───x───x───┘
           0.0     0.25    0.5     0.75    1.0

                            a

If something other than a unit square is desired, `lower_bound` and `upper_bound`
can be specified. Also, bases different than the default base 2 can be used,
e.g. `base=3` for a ternary grid.

In addition to the plain methods, there is a convenience layer for conversion
from the original coordinates
```julia-repl
julia> origcoord_to_grididx(grid; a=0.5, b=0.25)
(5, 2)

julia> origcoord_to_quantics(grid; a=0.5, b=0.25)
3-element Vector{Int64}:
 4
 1
 1
```
and also from grid indices
```julia-repl
julia> grididx_to_origcoord(grid; a=5, b=2)
(0.5, 0.25)

julia> grididx_to_quantics(grid; a=5, b=2)
3-element Vector{Int64}:
 4
 1
 1
```

For a simpler grid, we can just supply the resolution in each dimension:
```julia-repl
julia> boring_grid = NewDiscretizedGrid((3, 9))
NewDiscretizedGrid{2} with 8×512 = 4096 grid points
├─ Resolutions: (1: 3, 2: 9)
├─ Domain: unit square [0, 1)²
├─ Grid spacing: (Δ1 = 0.125, Δ2 = 0.001953125)
└─ Tensor train: 12 sites (uniform dimension 2)
```
In this case, variable names are automatically generated as `1`, `2`, etc.
"""
struct NewDiscretizedGrid{D}
    Rs::NTuple{D,Int}
    lower_bound::NTuple{D,Float64}
    upper_bound::NTuple{D,Float64}
    variablenames::NTuple{D,Symbol}
    base::Int
    indextable::Vector{Vector{Tuple{Symbol,Int}}}
    # Lookup table: lookup_table[variablename_index][bitnumber] -> (site_index, position_in_site)
    lookup_table::NTuple{D,Vector{Tuple{Int,Int}}}

    function NewDiscretizedGrid{D}(
        Rs, lower_bound, upper_bound, variablenames, base, indextable, includeendpoint
    ) where {D}
        @assert all(>=(0), Rs)
        @assert all(lower_bound .< upper_bound)
        @assert base > 1
        @assert all(R -> rangecheck_R(R; base), Rs)

        lookup_table = _build_lookup_table(Rs, indextable, variablenames, Val(D))

        includeendpoint = _to_tuple(Val(D), includeendpoint)
        lower_bound = _to_tuple(Val(D), lower_bound)
        upper_bound = _to_tuple(Val(D), upper_bound)

        upper_bound = _adjust_upper_bounds(
            upper_bound, lower_bound, includeendpoint, base, Rs, Val(D)
        )

        return new{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable, lookup_table)
    end
end

# Helper functions for constructor
function _build_lookup_table(Rs, indextable, variablenames, ::Val{D}) where D
    lookup_table = ntuple(D) do d
        Vector{Tuple{Int,Int}}(undef, Rs[d])
    end

    for (site_idx, quanticsindices) in pairs(indextable)
        for (pos_in_site, qindex) in pairs(quanticsindices)
            variablename, bitnumber = qindex
            var_idx = findfirst(==(variablename), variablenames)
            isnothing(var_idx) && continue
            lookup_table[var_idx][bitnumber] = (site_idx, pos_in_site)
        end
    end

    return lookup_table
end

function _adjust_upper_bounds(upper_bound, lower_bound, includeendpoint, base, Rs, ::Val{D}) where D
    return ntuple(D) do d
        if includeendpoint[d]
            upper_bound[d] + (upper_bound[d] - lower_bound[d]) / (base^Rs[d] - 1)
        else
            upper_bound[d]
        end
    end
end

function rangecheck_R(R::Int; base::Int=2)::Bool
    # For all gridindices from 1 to base^R to fit into an Int64,
    # we must have base ^ R <= typemax(Int)
    result = 1
    for _ in 1:R
        result <= typemax(Int) ÷ base || return false
        result *= base
    end
    return true
end

function NewDiscretizedGrid(
    Rs::NTuple{D,Int};
    lower_bound=default_lower_bound(Val(D)),
    upper_bound=default_upper_bound(Val(D)),
    base::Int=2,
    unfoldingscheme::Symbol=:interleaved,
    includeendpoint::Bool=false
) where {D}
    variablenames = ntuple(Symbol, D)
    indextable = _build_indextable(Rs, unfoldingscheme)

    return NewDiscretizedGrid{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable, includeendpoint)
end

function _build_indextable(Rs::NTuple{D,Int}, unfoldingscheme::Symbol) where D
    indextable = Vector{Tuple{Symbol,Int}}[]

    for bitnumber in 1:maximum(Rs)
        if unfoldingscheme === :interleaved
            _add_interleaved_indices!(indextable, Rs, bitnumber)
        elseif unfoldingscheme === :fused
            _add_fused_indices!(indextable, Rs, bitnumber)
        else
            throw(ArgumentError(lazy"""Unfolding scheme $unfoldingscheme not supported. Use :interleaved or :fused.
            If you need a different scheme, please use the NewDiscretizedGrid(variablenames::NTuple{D,Symbol}, indextable::Vector{Vector{Tuple{Symbol,Int}}}) constructor."""))
        end
    end

    return indextable
end

function _add_interleaved_indices!(indextable, Rs::NTuple{D,Int}, bitnumber) where D
    for d in 1:D
        bitnumber ∈ 1:Rs[d] || continue
        variablename = Symbol(d)
        qindex = (variablename, bitnumber)
        push!(indextable, [qindex])
    end
end

function _add_fused_indices!(indextable, Rs::NTuple{D,Int}, bitnumber) where D
    indices_bitnumber = Tuple{Symbol,Int}[]
    for d in 1:D
        bitnumber ∈ 1:Rs[d] || continue
        variablename = Symbol(d)
        qindex = (variablename, bitnumber)
        push!(indices_bitnumber, qindex)
    end
    if !isempty(indices_bitnumber)
        push!(indextable, indices_bitnumber)
    end
end

function NewDiscretizedGrid{D}(
    R::Int,
    lower_bound=default_lower_bound(Val(D)),
    upper_bound=default_upper_bound(Val(D));
    kwargs...
) where {D}
    return NewDiscretizedGrid(ntuple(Returns(R), D); lower_bound, upper_bound, kwargs...)
end

function NewDiscretizedGrid(
    variablenames::NTuple{D,Symbol},
    indextable::Vector{Vector{Tuple{Symbol,Int}}};
    lower_bound=default_lower_bound(Val(D)),
    upper_bound=default_upper_bound(Val(D)),
    base::Int=2,
    includeendpoint::Bool=false
) where D
    @assert all(Iterators.flatten(indextable)) do index
        first(index) ∈ variablenames
    end

    Rs = Tuple(map(variablenames) do variablename
        count(index -> first(index) == variablename, Iterators.flatten(indextable))
    end)

    return NewDiscretizedGrid{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable, includeendpoint)
end

default_lower_bound(::Val{D}) where D = ntuple(Returns(0.0), D)

default_upper_bound(::Val{D}) where D = ntuple(Returns(1.0), D)

Base.ndims(::NewDiscretizedGrid{D}) where D = D

Base.length(g::NewDiscretizedGrid) = length(g.indextable)

function sitedim(g::NewDiscretizedGrid, site::Int)::Int
    @assert site ∈ eachindex(g.indextable)
    return g.base^length(g.indextable[site])
end

function quantics_to_grididx(g::NewDiscretizedGrid{D}, quantics::AbstractVector{Int}) where D
    @assert length(quantics) == length(g)
    @assert all(site -> quantics[site] ∈ 1:sitedim(g, site), eachindex(quantics))

    result = if g.base == 2
        _quantics_to_grididx_base2(g, quantics)
    else
        _quantics_to_grididx_general(g, quantics)
    end

    return _convert_to_scalar_if_possible(result)
end

function _quantics_to_grididx_general(g::NewDiscretizedGrid{D}, quantics) where D
    base = g.base

    return ntuple(D) do d
        grididx = 1
        R_d = g.Rs[d]

        for bitnumber in 1:R_d
            site_idx, pos_in_site = g.lookup_table[d][bitnumber]
            quantics_val = quantics[site_idx]
            site_len = length(g.indextable[site_idx])

            temp = quantics_val - 1
            for _ in 1:(site_len-pos_in_site)
                temp = div(temp, base)
            end
            digit = temp % base

            grididx += digit * base^(R_d - bitnumber)
        end
        grididx
    end
end

function _quantics_to_grididx_base2(g::NewDiscretizedGrid{D}, quantics) where D
    return ntuple(D) do d
        grididx = 0
        R_d = g.Rs[d]

        for bitnumber in 1:R_d
            site_idx, pos_in_site = g.lookup_table[d][bitnumber]
            bit_position = length(g.indextable[site_idx]) - pos_in_site
            digit = ((quantics[site_idx] - 1) >> bit_position) & 1
            grididx |= digit << (R_d - bitnumber)
        end
        grididx + 1
    end
end

function grididx_to_quantics(g::NewDiscretizedGrid{D}, grididx) where D
    grididx_tuple = _to_tuple(Val(D), grididx)

    @assert length(grididx_tuple) == D lazy"Grid index must have dimension $D, got $(length(grididx_tuple))"
    @assert all(1 ≤ grididx_tuple[d] ≤ g.base^g.Rs[d] for d in 1:D) "Grid index out of bounds"

    result = ones(Int, length(g.indextable))

    if g.base == 2
        _grididx_to_quantics_base2!(result, g, grididx_tuple)
    else
        _grididx_to_quantics_general!(result, g, grididx_tuple)
    end

    return result
end

function _grididx_to_quantics_general!(result::Vector{Int}, g::NewDiscretizedGrid{D}, grididx::NTuple{D,Int}) where D
    base = g.base

    @inbounds for d in 1:D
        zero_based_idx = grididx[d] - 1
        R_d = g.Rs[d]

        for bitnumber in 1:R_d
            site_idx, pos_in_site = g.lookup_table[d][bitnumber]
            site_length = length(g.indextable[site_idx])

            bit_position = R_d - bitnumber
            digit = (zero_based_idx ÷ (base^bit_position)) % base

            power = site_length - pos_in_site
            result[site_idx] += digit * (base^power)
        end
    end
end

function _grididx_to_quantics_base2!(result::Vector{Int}, g::NewDiscretizedGrid{D}, grididx::NTuple{D,Int}) where D
    @inbounds for d in 1:D
        zero_based_idx = grididx[d] - 1
        R_d = g.Rs[d]

        for bitnumber in 1:R_d
            site_idx, pos_in_site = g.lookup_table[d][bitnumber]
            site_length = length(g.indextable[site_idx])

            bit_position = R_d - bitnumber
            digit = (zero_based_idx >> bit_position) & 1

            power = site_length - pos_in_site
            result[site_idx] += digit << power
        end
    end
end

grid_min(g::NewDiscretizedGrid) = _convert_to_scalar_if_possible(g.lower_bound)

grid_max(g::NewDiscretizedGrid) = _convert_to_scalar_if_possible(g.upper_bound .- grid_step(g))

grid_step(g::NewDiscretizedGrid) = _convert_to_scalar_if_possible(
    (upper_bound(g) .- lower_bound(g)) ./ (g.base .^ g.Rs),
)

upper_bound(g::NewDiscretizedGrid) = _convert_to_scalar_if_possible(g.upper_bound)

lower_bound(g::NewDiscretizedGrid) = _convert_to_scalar_if_possible(g.lower_bound)

function grididx_to_origcoord(g::NewDiscretizedGrid{D}, index) where {D}
    index_tuple = _to_tuple(Val(D), index)
    @assert all(1 .<= index .<= (g.base .^ g.Rs)) lazy"Grid-index $index out of bounds [1, $(g.base .^ g.Rs)]"
    return _convert_to_scalar_if_possible((index .- 1) .* grid_step(g) .+ grid_min(g))
end

function origcoord_to_grididx(g::NewDiscretizedGrid{D}, coordinate) where {D}
    coord_tuple = _to_tuple(Val(D), coordinate)

    bounds_lower = lower_bound(g)
    bounds_upper = upper_bound(g)
    @assert all(bounds_lower .<= coord_tuple .<= bounds_upper) "Coordinate $coord_tuple out of bounds [$bounds_lower, $bounds_upper]"

    min_vals = grid_min(g)
    step_vals = grid_step(g)
    raw_indices = @. round(Int, (coord_tuple - min_vals) / step_vals + 1)
    min_indices = ntuple(Returns(1), D)
    max_indices = g.base .^ g.Rs
    clamped_indices = clamp.(raw_indices, min_indices, max_indices)

    return _convert_to_scalar_if_possible(clamped_indices)
end

function origcoord_to_quantics(g::NewDiscretizedGrid{D}, coordinate) where {D}
    coord_tuple = _to_tuple(Val(D), coordinate)
    grid_idx = origcoord_to_grididx(g, coord_tuple)
    return grididx_to_quantics(g, grid_idx)
end

function quantics_to_origcoord(g::NewDiscretizedGrid{D}, quantics::AbstractVector{Int}) where {D}
    grid_idx = quantics_to_grididx(g, quantics)
    return grididx_to_origcoord(g, grid_idx)
end

function localdimensions(g::NewDiscretizedGrid)::Vector{Int}
    return g.base .^ length.(g.indextable)
end

function quanticsfunction(::Type{T}, g::NewDiscretizedGrid, f::F)::Function where {T,F<:Function}
    function wrapped_function(quantics)::T
        coords = quantics_to_origcoord(g, quantics)
        if coords isa Tuple
            return f(coords...)
        else
            return f(coords)
        end
    end
    return wrapped_function
end
grid_origin(g::NewDiscretizedGrid) = lower_bound(g)

# Keyword argument convenience functions

function _handle_kwargs_input(g::NewDiscretizedGrid{D}; kwargs...) where {D}
    provided_keys = keys(kwargs)
    expected_keys = g.variablenames
    @assert Set(provided_keys) == Set(expected_keys) "Expected keyword arguments $(expected_keys), got $(tuple(provided_keys...))"

    @assert all(v -> v isa Real, values(kwargs)) "All keyword argument values must be Real numbers"

    return ntuple(D) do d
        variablename = g.variablenames[d]
        kwargs[variablename]
    end
end

function origcoord_to_grididx(g::NewDiscretizedGrid; kwargs...)
    coordinate = _handle_kwargs_input(g; kwargs...)
    return origcoord_to_grididx(g, coordinate)
end

function origcoord_to_quantics(g::NewDiscretizedGrid; kwargs...)
    coordinate = _handle_kwargs_input(g; kwargs...)
    return origcoord_to_quantics(g, coordinate)
end

function grididx_to_origcoord(g::NewDiscretizedGrid; kwargs...)
    index = _handle_kwargs_input(g; kwargs...)
    return grididx_to_origcoord(g, index)
end

function grididx_to_quantics(g::NewDiscretizedGrid; kwargs...)
    index = _handle_kwargs_input(g; kwargs...)
    return grididx_to_quantics(g, index)
end

function Base.show(io::IO, ::MIME"text/plain", g::NewDiscretizedGrid{D}) where D
    print(io, "NewDiscretizedGrid{$D}")

    # Grid resolution and total points
    total_points = prod(g.base .^ g.Rs)
    if D == 1
        print(io, " with $(g.base^g.Rs[1]) grid points")
    else
        print(io, " with $(join(g.base .^ g.Rs, "×")) = $total_points grid points")
    end

    # Variable names (if meaningful)
    if any(name -> !startswith(string(name), r"^\d+$"), g.variablenames)
        var_str = join(g.variablenames, ", ")
        print(io, "\n├─ Variables: ($var_str)")
    end

    # Resolution per dimension
    if D == 1
        print(io, "\n├─ Resolution: $(g.Rs[1]) bits")
    else
        res_str = join(["$(g.variablenames[i]): $(g.Rs[i])" for i in 1:D], ", ")
        print(io, "\n├─ Resolutions: ($res_str)")
    end

    # Bounds (only show if not default unit interval/square/cube)
    default_lower = ntuple(Returns(0.0), D)
    default_upper = ntuple(Returns(1.0), D)
    if g.lower_bound != default_lower || any(abs.(g.upper_bound .- default_upper) .> 1e-10)
        if D == 1
            print(io, "\n├─ Domain: [$(g.lower_bound[1]), $(g.upper_bound[1]))")
        else
            bounds_str = join(["[$(g.lower_bound[i]), $(g.upper_bound[i]))" for i in 1:D], " × ")
            print(io, "\n├─ Domain: $bounds_str")
        end

        # Grid spacing
        step_vals = grid_step(g)
        if D == 1
            print(io, "\n├─ Grid spacing: $(step_vals)")
        else
            step_str = join(["Δ$(g.variablenames[i]) = $(step_vals[i])" for i in 1:D], ", ")
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
            step_str = join(["Δ$(g.variablenames[i]) = $(step_vals[i])" for i in 1:D], ", ")
            print(io, "\n├─ Grid spacing: ($step_str)")
        end
    end

    # Base (only show if not binary)
    if g.base != 2
        print(io, "\n├─ Base: $(g.base)")
    end

    # Tensor structure summary
    num_sites = length(g.indextable)
    site_dims = [g.base^length(site) for site in g.indextable]
    max_bond_dim = maximum(site_dims)

    print(io, "\n└─ Tensor train: $num_sites sites")
    if all(d -> d == site_dims[1], site_dims)
        print(io, " (uniform dimension $(site_dims[1]))")
    else
        print(io, " (dimensions: $(join(site_dims, "-")))")
    end
end

function Base.show(io::IO, g::NewDiscretizedGrid{D}) where D
    total_points = prod(g.base .^ g.Rs)
    if D == 1
        print(io, "NewDiscretizedGrid{$D}($(g.base^g.Rs[1]) points)")
    else
        print(io, "NewDiscretizedGrid{$D}($(join(g.base .^ g.Rs, "×")) points)")
    end
end
