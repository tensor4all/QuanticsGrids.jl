"""
    DiscretizedGrid{D}

A discretized grid structure for D-dimensional grids with variable resolution,
supporting efficient conversion between quantics, grid indices, and original coordinates.
A `DiscretizedGrid` instance is intended to undergird a quantics tensor train
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
julia> grid = DiscretizedGrid((:a, :b), [[(:a, 1), (:b, 2)], [(:a, 2)], [(:b, 1), (:a, 3)]])
DiscretizedGrid{2} with 8×4 = 32 grid points
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
e.g. `base=3` for a ternary grid, or `base=(2, 6)` for per-variable bases. When
passing a tuple, the base order matches `variablenames`.

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
julia> boring_grid = DiscretizedGrid((3, 9))
DiscretizedGrid{2} with 8×512 = 4096 grid points
├─ Resolutions: (1: 3, 2: 9)
├─ Domain: unit square [0, 1)²
├─ Grid spacing: (Δ1 = 0.125, Δ2 = 0.001953125)
└─ Tensor train: 12 sites (uniform dimension 2)
```
In this case, variable names are automatically generated as `1`, `2`, etc.
"""
struct DiscretizedGrid{D} <: Grid{D}
    discretegrid::InherentDiscreteGrid{D}
    lower_bound::NTuple{D,Float64}
    upper_bound::NTuple{D,Float64}
    includeendpoint::NTuple{D,Bool}

    function DiscretizedGrid{D}(
        Rs, lower_bound, upper_bound, variablenames, base, indextable, includeendpoint
    ) where {D}
        lower_bound = _to_tuple(Val(D), lower_bound)
        upper_bound = _to_tuple(Val(D), upper_bound)
        base = _to_tuple(Val(D), base)
        includeendpoint = _to_tuple(Val(D), includeendpoint)
        for d in 1:D
            if !(lower_bound[d] < upper_bound[d])
                throw(ArgumentError(lazy"Got (lower_bound[$d], upper_bound[$d]) = $((lower_bound[d], upper_bound[d])). Each lower bound needs to be strictly less than the corresponding upper bound."))
            end
        end

        discretegrid = InherentDiscreteGrid(variablenames, indextable; base)

        upper_bound = _adjust_upper_bounds(
            upper_bound, lower_bound, includeendpoint, base, Rs, Val(D)
        )

        return new{D}(discretegrid, lower_bound, upper_bound, includeendpoint)
    end
end

# ============================================================================
# Helper/utility functions
# ============================================================================

function _check_bounds_dim(::Val{D}, lower_bound, upper_bound) where {D}
    if lower_bound isa NTuple && length(lower_bound) != D
        throw(ArgumentError(lazy"Got lower_bound with length $(length(lower_bound)); expected $D for DiscretizedGrid{$D}."))
    end
    if upper_bound isa NTuple && length(upper_bound) != D
        throw(ArgumentError(lazy"Got upper_bound with length $(length(upper_bound)); expected $D for DiscretizedGrid{$D}."))
    end
end

function _adjust_upper_bounds(upper_bound, lower_bound, includeendpoint, base, Rs, ::Val{D}) where D
    base = _to_tuple(Val(D), base)
    includeendpoint = _to_tuple(Val(D), includeendpoint)
    for d in 1:D
        if iszero(Rs[d]) && includeendpoint[d]
            throw(ArgumentError(lazy"Got Rs[$d] = 0 and includeendpoint[$d] = true. This is not allowed."))
        end
    end

    return ntuple(D) do d
        if includeendpoint[d]
            base_d = base[d]
            upper_bound[d] + (upper_bound[d] - lower_bound[d]) / (base_d^Rs[d] - 1)
        else
            upper_bound[d]
        end
    end
end

default_lower_bound(::Val{D}) where D = _to_tuple(Val(D), 0.0)

default_upper_bound(::Val{D}) where D = _to_tuple(Val(D), 1.0)

function _handle_kwargs_input(g::DiscretizedGrid{D}; kwargs...) where {D}
    provided_keys = keys(kwargs)
    expected_keys = grid_variablenames(g)
    if !(Set(provided_keys) == Set(expected_keys))
        throw(ArgumentError(lazy"Expected keyword arguments $(expected_keys), got $(tuple(provided_keys...))"))
    end

    if !all(v -> v isa Real, values(kwargs))
        throw(ArgumentError(lazy"Got kwargs = $kwargs. All keyword argument values must be Real numbers"))
    end

    return ntuple(D) do d
        variablename = grid_variablenames(g)[d]
        kwargs[variablename]
    end
end

# ============================================================================
# Constructors
# ============================================================================

function DiscretizedGrid(
    variablenames::NTuple{D,Symbol},
    Rs::NTuple{D,Int};
    lower_bound=default_lower_bound(Val(D)),
    upper_bound=default_upper_bound(Val(D)),
    base=2,
    unfoldingscheme::Symbol=:fused,
    includeendpoint=false
) where {D}
    indextable = _build_indextable(variablenames, Rs, unfoldingscheme)

    return DiscretizedGrid{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable, includeendpoint)
end

function DiscretizedGrid(Rs::NTuple{D,Int}; variablenames=ntuple(Symbol, D), kwargs...) where {D}
    return DiscretizedGrid(variablenames, Rs; kwargs...)
end

function DiscretizedGrid{D}(
    R,
    lower_bound=default_lower_bound(Val(D)),
    upper_bound=default_upper_bound(Val(D));
    kwargs...
) where {D}
    _check_bounds_dim(Val(D), lower_bound, upper_bound)
    return DiscretizedGrid(_to_tuple(Val(D), R); lower_bound, upper_bound, kwargs...)
end

function DiscretizedGrid(
    R,
    lower_bound::NTuple{D,Real},
    upper_bound::NTuple{D,Real};
    kwargs...
) where {D}
    _check_bounds_dim(Val(D), lower_bound, upper_bound)
    return DiscretizedGrid(_to_tuple(Val(D), R); lower_bound, upper_bound, kwargs...)
end

function DiscretizedGrid(
    variablenames::NTuple{D,Symbol},
    indextable::Vector{Vector{Tuple{Symbol,Int}}};
    lower_bound=default_lower_bound(Val(D)),
    upper_bound=default_upper_bound(Val(D)),
    base=2,
    includeendpoint=false
) where D
    Rs = map(variablenames) do variablename
        count(index -> first(index) == variablename, Iterators.flatten(indextable))
    end

    return DiscretizedGrid{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable, includeendpoint)
end

function DiscretizedGrid(R::Int, lower_bound::Real, upper_bound::Real; kwargs...)
    return DiscretizedGrid{1}(R, lower_bound, upper_bound; kwargs...)
end

# ============================================================================
# Basic property accessor functions
# ============================================================================

Base.ndims(::DiscretizedGrid{D}) where D = D

Base.length(g::DiscretizedGrid) = length(grid_indextable(g))

grid_Rs(g::DiscretizedGrid{D}) where D = grid_Rs(g.discretegrid)

grid_indextable(g::DiscretizedGrid{D}) where D = grid_indextable(g.discretegrid)

grid_base(g::DiscretizedGrid{D}) where D = grid_base(g.discretegrid)

grid_bases(g::DiscretizedGrid{D}) where D = grid_bases(g.discretegrid)

grid_variablenames(g::DiscretizedGrid{D}) where D = grid_variablenames(g.discretegrid)

upper_bound(g::DiscretizedGrid) = _convert_to_scalar_if_possible(g.upper_bound)

lower_bound(g::DiscretizedGrid) = _convert_to_scalar_if_possible(g.lower_bound)

grid_origin(g::DiscretizedGrid) = lower_bound(g)

function sitedim(g::DiscretizedGrid, site::Int)::Int
    if !(site ∈ eachindex(grid_indextable(g)))
        throw(DomainError(site, lazy"Site index out of bounds [1, $(length(grid_indextable(g)))]"))
    end
    return sitedim(g.discretegrid, site)
end

# ============================================================================
# Grid coordinate functions
# ============================================================================

grid_min(g::DiscretizedGrid) = _convert_to_scalar_if_possible(g.lower_bound)

grid_max(g::DiscretizedGrid) = _convert_to_scalar_if_possible(g.upper_bound .- grid_step(g))

grid_step(g::DiscretizedGrid) = _convert_to_scalar_if_possible(
    (upper_bound(g) .- lower_bound(g)) ./ (grid_bases(g) .^ grid_Rs(g)),
)

function grid_origcoords(g::DiscretizedGrid, d::Int)
    if !(1 ≤ d ≤ ndims(g))
        throw(DomainError(d, lazy"Dimension $d out of bounds [1, $(ndims(g))]."))
    end
    start = grid_min(g)[d]
    stop = grid_max(g)[d]
    base_d = grid_bases(g)[d]
    length = base_d^grid_Rs(g)[d]
    return range(start, stop, length)
end

function grid_origcoords(g::DiscretizedGrid, variablename::Symbol)
    d = findfirst(==(variablename), grid_variablenames(g))
    isnothing(d) && throw(ArgumentError(lazy"Variable name :$variablename not found in grid. Available variables: $(grid_variablenames(g))"))
    return grid_origcoords(g, d)
end

# ============================================================================
# Core conversion functions
# ============================================================================

function quantics_to_grididx(g::DiscretizedGrid, quantics::AbstractVector{Int})
    return quantics_to_grididx(g.discretegrid, quantics)
end

function grididx_to_quantics(g::DiscretizedGrid, grididx)
    return grididx_to_quantics(g.discretegrid, grididx)
end

function grididx_to_origcoord(g::DiscretizedGrid{D}, index) where {D}
    index_tuple = _to_tuple(Val(D), index)
    for d in 1:D
        base_d = grid_bases(g)[d]
        if !(1 <= index_tuple[d] <= base_d^grid_Rs(g)[d])
            throw(DomainError(index_tuple[d], lazy"Grid index out of bounds [1, $(base_d ^ grid_Rs(g)[d])]."))
        end
    end

    res = ntuple(D) do d
        base_d = grid_bases(g)[d]
        step_d = (g.upper_bound[d] - g.lower_bound[d]) / (base_d^grid_Rs(g)[d])
        g.lower_bound[d] + (index_tuple[d] - 1) * step_d
    end

    return _convert_to_scalar_if_possible(res)
end

function origcoord_to_grididx(g::DiscretizedGrid{D}, coordinate) where {D}
    coord_tuple = _to_tuple(Val(D), coordinate)

    bounds_lower = lower_bound(g)
    bounds_upper = upper_bound(g)
    for d in 1:D
        if !(bounds_lower[d] <= coord_tuple[d] <= bounds_upper[d])
            throw(DomainError(coord_tuple[d], lazy"Coordinate out of bounds [$(bounds_lower[d]), $(bounds_upper[d])] (dimension $d)."))
        end
    end

    steps = grid_step(g)
    indices = ntuple(D) do d
        target = coord_tuple[d]
        step_d = steps[d]

        continuous_idx = (target - bounds_lower[d]) / step_d + 1

        discrete_idx = round(Int, continuous_idx)
        base_d = grid_bases(g)[d]
        clamp(discrete_idx, 1, base_d^grid_Rs(g)[d])
    end

    return _convert_to_scalar_if_possible(indices)
end

function origcoord_to_quantics(g::DiscretizedGrid{D}, coordinate) where {D}
    coord_tuple = _to_tuple(Val(D), coordinate)
    grid_idx = origcoord_to_grididx(g, coord_tuple)
    return grididx_to_quantics(g, grid_idx)
end

function quantics_to_origcoord(g::DiscretizedGrid{D}, quantics::AbstractVector{Int}) where {D}
    grid_idx = quantics_to_grididx(g, quantics)
    return grididx_to_origcoord(g, grid_idx)
end

# ============================================================================
# Keyword argument convenience functions
# ============================================================================

function origcoord_to_grididx(g::DiscretizedGrid; kwargs...)
    coordinate = _handle_kwargs_input(g; kwargs...)
    return origcoord_to_grididx(g, coordinate)
end

function origcoord_to_quantics(g::DiscretizedGrid; kwargs...)
    coordinate = _handle_kwargs_input(g; kwargs...)
    return origcoord_to_quantics(g, coordinate)
end

function grididx_to_origcoord(g::DiscretizedGrid; kwargs...)
    index = _handle_kwargs_input(g; kwargs...)
    return grididx_to_origcoord(g, index)
end

function grididx_to_quantics(g::DiscretizedGrid; kwargs...)
    index = _handle_kwargs_input(g; kwargs...)
    return grididx_to_quantics(g, index)
end

# ============================================================================
# Other utility functions
# ============================================================================

function localdimensions(g::DiscretizedGrid)::Vector{Int}
    return copy(g.discretegrid.sitedims)
end

function quanticsfunction(::Type{T}, g::Grid, f::F)::Function where {T,F<:Function}
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

# ============================================================================
# Display/show methods
# ============================================================================

function Base.show(io::IO, ::MIME"text/plain", g::DiscretizedGrid{D}) where D
    print(io, "DiscretizedGrid{$D}")

    # Grid resolution and total points
    points_per_dim = big.(grid_bases(g)) .^ grid_Rs(g)
    total_points = prod(points_per_dim)
    if D <= 1
        print(io, " with $total_points grid point" * (isone(total_points) ? "" : "s"))
    else
        print(io, " with $(join(points_per_dim, " × ")) = $total_points grid points")
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
            upper_display = g.includeendpoint[1] ? grid_max(g) : g.upper_bound[1]
            upper_bracket = g.includeendpoint[1] ? "]" : ")"
            print(io, "\n├─ Domain: [$(g.lower_bound[1]), $(upper_display)$(upper_bracket)")
        else
            upper_display = ntuple(D) do i
                g.includeendpoint[i] ? grid_max(g)[i] : g.upper_bound[i]
            end
            bounds_str = join([
                    "[$(g.lower_bound[i]), $(upper_display[i])$(g.includeendpoint[i] ? "]" : ")")"
                    for i in 1:D
                ], " × ")
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
    sitedims = Int[sitedim(g, site) for site in 1:num_sites]

    print(io, "\n└─ Tensor train: $num_sites sites")
    if !isempty(sitedims)
        if allequal(sitedims)
            print(io, " (uniform dimension $(sitedims[1]))")
        else
            print(io, " (dimensions: $(join(sitedims, "-")))")
        end
    end
end

function Base.show(io::IO, g::DiscretizedGrid{D}) where D
    points_per_dim = grid_bases(g) .^ grid_Rs(g)
    if D == 1
        print(io, "DiscretizedGrid{$D}($(points_per_dim[1]) points)")
    else
        print(io, "DiscretizedGrid{$D}($(join(points_per_dim, "×")) points)")
    end
end
