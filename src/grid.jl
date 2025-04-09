
# Grid for d-dimensional space
abstract type Grid{d} end


_convert_to_scalar_if_possible(x) = x
_convert_to_scalar_if_possible(x::NTuple{1,T}) where {T} = first(x)

_to_tuple(::Val{d}, x::NTuple{d,T}) where {d,T} = x
_to_tuple(::Val{d}, x) where {d} = ntuple(i -> x, d)

function digitmax(d, base, unfoldingscheme::Symbol)
    if unfoldingscheme === :fused
        return base^d
    else
        return base
    end
end

function quanticslength(R, d, unfoldingscheme::Symbol)
    if unfoldingscheme === :fused
        return R
    else
        return R * d
    end
end

function localdimensions(g::Grid{d}) where {d}
    return fill(
        digitmax(d, g.base, g.unfoldingscheme),
        quanticslength(g.R, d, g.unfoldingscheme)
    )
end

function _rangecheck_R(R; base=2)::Int
    base * (BigInt(base)^R - 1) ÷ (base - 1) <= typemax(Int) ||
        error("R too large for base $base")
end

#===
Conversion between grid index, original coordinate, quantics

* grid index => original coordinate (generic)
* grid index => quantics (generic)

* original coordinate => grid index (implemented for each grid type)
* original coordinate => quantics (generic)

* quantics => grididx (XXXX)
* quantics => origcoord (generic)
===#

"""
grid index => original coordinate
"""
function grididx_to_origcoord(g::Grid{d}, index::NTuple{d,Int}) where {d}
    all(1 .<= index .<= g.base^g.R) || error("1 <= {index} <= g.base^g.R")
    return _convert_to_scalar_if_possible((index .- 1) .* grid_step(g) .+ grid_min(g))
end

"""
grid index => original coordinate
"""
grididx_to_origcoord(g::Grid{1}, index::Int) = grididx_to_origcoord(g, (index,))

"""
grid index => quantics
"""
function grididx_to_quantics(g::Grid{d}, grididx::NTuple{d,Int}) where {d}
    if g.unfoldingscheme === :fused
        return index_to_quantics_fused(grididx, numdigits=g.R, base=g.base)
    else
        return fused_to_interleaved(
            index_to_quantics_fused(grididx, numdigits=g.R, base=g.base),
            d,
            base=g.base,
        )
    end
end

"""
grid index => quantics
"""
grididx_to_quantics(g::Grid{1}, grididx::Int) = grididx_to_quantics(g, (grididx,))


"""
original coordinate => quantics
"""
function origcoord_to_quantics(g::Grid{d}, coordinate) where {d}
    return grididx_to_quantics(g, origcoord_to_grididx(g, _to_tuple(Val(d), coordinate)))
end


"""
quantics => grid index
"""
function quantics_to_grididx(g::Grid{d}, bitlist) where {d}
    (all(1 .<= bitlist) && all(bitlist .<= digitmax(d, g.base, g.unfoldingscheme))) ||
        error("Range error for bitlist: $bitlist")
    quanticslength(g.R, d, g.unfoldingscheme) == length(bitlist) ||
        error("Length error for bitlist: $bitlist")

    return _convert_to_scalar_if_possible(
        quantics_to_index(
            bitlist;
            base=g.base,
            dims=Val(d),
            unfoldingscheme=g.unfoldingscheme,
        ),
    )
end


"""
quantics => original coordinate system
"""
function quantics_to_origcoord(g::Grid{d}, bitlist) where {d}
    return _convert_to_scalar_if_possible(
        grididx_to_origcoord(g, quantics_to_grididx(g, bitlist)),
    )
end


"""
Make a wrapper function that takes a bitlist as input
"""
function quanticsfunction(::Type{T}, g::Grid{d}, f::Function)::Function where {T,d}
    function _f(bitlist)::T
        return f(quantics_to_origcoord(g, bitlist)...)
    end
    return _f
end


@doc raw"""
The InherentDiscreteGrid struct represents a grid for inherently discrete data.
The grid contains values at specific,
equally spaced points, but these values do not represent discretized versions
of continuous data. Instead, they represent individual data points that are
inherently discrete.
The linear size of the mesh is ``base^R``, where ``base`` defaults to 2.
"""
struct InherentDiscreteGrid{d} <: Grid{d}
    R::Int
    origin::NTuple{d,Int}
    base::Int
    unfoldingscheme::Symbol
    step::NTuple{d,Int}

    function InherentDiscreteGrid{d}(
        R::Int,
        origin::Union{NTuple{d,Int},Int};
        base::Integer=2,
        unfoldingscheme::Symbol=:fused,
        step::Union{NTuple{d,Int},Int}=1,
    ) where {d}
        _rangecheck_R(R; base=base)
        unfoldingscheme in (:fused, :interleaved) ||
            error("Invalid unfolding scheme: $unfoldingscheme")
        origin_ = origin isa Int ? ntuple(i -> origin, d) : origin
        step_ = step isa Int ? ntuple(i -> step, d) : step
        new(R, origin_, base, unfoldingscheme, step_)
    end

end

grid_min(grid::InherentDiscreteGrid) = _convert_to_scalar_if_possible(grid.origin)
grid_max(grid::InherentDiscreteGrid) = _convert_to_scalar_if_possible(grid.origin .+ grid_step(grid) .* (grid.base^grid.R .- 1))
grid_step(grid::InherentDiscreteGrid{d}) where {d} =
    _convert_to_scalar_if_possible(grid.step)

"""
Create a grid for inherently discrete data with origin at 1
"""
function InherentDiscreteGrid{d}(
    R::Int;
    base::Integer=2,
    step::Union{NTuple{d,Int},Int}=1,
    unfoldingscheme::Symbol=:fused,
) where {d}
    InherentDiscreteGrid{d}(
        R,
        1;
        base=base,
        unfoldingscheme=unfoldingscheme,
        step=step,
    )
end


"""
Create a grid for inherently discrete data with origin at an arbitrary point
"""
#function InherentDiscreteGrid{d}(
#R::Int, origin::Union{NTuple{d,Int},Int};
#base=2, unfoldingscheme::UnfoldingSchemes.UnfoldingScheme=UnfoldingSchemes.fused,
#) where {d}
#InherentDiscreteGrid{d}(
#R, origin; base=base, unfoldingscheme=unfoldingscheme, step=step
#)
#end


"""
Convert a coordinate in the original coordinate system to the corresponding grid index
"""
function origcoord_to_grididx(
    g::InherentDiscreteGrid,
    coordinate::Union{Int,NTuple{N,Int}},
) where {N}
    return _convert_to_scalar_if_possible(
        div.(coordinate .- grid_min(g), grid_step(g)) .+ 1,
    )
end


@doc raw"""
The DiscretizedGrid struct represents a grid for discretized continuous data.
This is used for data that is originally continuous,
but has been discretized for computational purposes.
The grid contains values at specific, equally spaced points, which represent discrete
approximations of the original continuous data.
"""
struct DiscretizedGrid{d} <: Grid{d}
    R::Int
    lower_bound::NTuple{d,Float64}
    upper_bound::NTuple{d,Float64}
    base::Int
    unfoldingscheme::Symbol
    includeendpoint::Bool

    function DiscretizedGrid{d}(
        R::Int,
        lower_bound,
        upper_bound;
        base::Integer=2,
        unfoldingscheme::Symbol=:fused,
        includeendpoint::Bool=false,
    ) where {d}
        _rangecheck_R(R; base=base)
        unfoldingscheme in (:fused, :interleaved) ||
            error("Invalid unfolding scheme: $unfoldingscheme")
        return new(
            R,
            _to_tuple(Val(d), lower_bound),
            _to_tuple(Val(d), upper_bound),
            base,
            unfoldingscheme,
            includeendpoint,
        )
    end
end

"""
    grid_min(g::DiscretizedGrid)

Returns the grid point with minimal coordinate values.
This is equivalent to [`lower_bound`](@ref).
The return value is scalar for 1D grids, and a `Tuple` otherwise.
"""
grid_min(g::DiscretizedGrid) = _convert_to_scalar_if_possible(g.lower_bound)

"""
    grid_max(g::DiscretizedGrid)

Returns the grid point with maximal coordinate values.
 - If `includeendpoint=false` during creation of the grid, this value is dependent on grid resolution.
 - If `includeendpoint=true` during creation of the grid, this function is equivalent to [`upper_bound`](@ref).
The return value is scalar for 1D grids, and a `Tuple` otherwise.
"""
grid_max(g::DiscretizedGrid) = g.includeendpoint ? _convert_to_scalar_if_possible(g.upper_bound) : _convert_to_scalar_if_possible(g.upper_bound .- grid_step(g))

"""
    grid_step(g::DiscretizedGrid)

Returns the distance between adjacent grid points in each dimension.
The return value is scalar for 1D grids, and a `Tuple` otherwise.
"""
grid_step(g::DiscretizedGrid{d}) where {d} = _convert_to_scalar_if_possible(
    g.includeendpoint ?
    (g.upper_bound .- g.lower_bound) ./ (g.base^g.R .- 1) :
    (g.upper_bound .- g.lower_bound) ./ (g.base^g.R)
)
"""
    upper_bound(g::DiscretizedGrid)

Returns the upper bound of the grid coordinates, as passed to the constructor during grid creation.
 - If `includeendpoint=false` during grid creation, this function returns a point that is one grid spacing beyond the last grid point (which can be obtained through [`grid_max`](@ref)).
 - If `includeendpoint=true` during grid creation, this function returns the point with maximal coordinate values. This is equivalent to [`grid_max`](@ref).
The return value is scalar for 1D grids, and a `Tuple` otherwise.
"""
upper_bound(g::DiscretizedGrid) = _convert_to_scalar_if_possible(g.upper_bound)

"""
    grid_min(g::DiscretizedGrid)

Returns the grid point with minimal coordinate values, as passed to the constructor during grid creation.
This is equivalent to [`grid_min`](@ref).
The return value is scalar for 1D grids, and a `Tuple` otherwise.
"""
lower_bound(g::DiscretizedGrid) = _convert_to_scalar_if_possible(g.lower_bound)

function DiscretizedGrid{d}(
    R::Int;
    base=2,
    unfoldingscheme::Symbol=:fused,
) where {d}
    return DiscretizedGrid{d}(
        R,
        ntuple(i -> 0.0, d),
        ntuple(i -> 1.0, d);
        base=base,
        unfoldingscheme=unfoldingscheme,
    )
end


"""
Convert a coordinate in the original coordinate system to the corresponding grid index
"""
function origcoord_to_grididx(g::DiscretizedGrid, coordinate::NTuple{N,Float64}) where {N}
    all(lower_bound(g) .<= coordinate .<= upper_bound(g)) ||
        error("Bound Error: $(coordinate), lower_bound=$(lower_bound(g)), upper_bound=$(upper_bound(g))")
    clip(x) = max.(1, min.(x, g.base^g.R))
    return _convert_to_scalar_if_possible(
        ((coordinate .- grid_min(g)) ./ grid_step(g) .+ 1) .|> round .|> Int,
    )
end

function origcoord_to_grididx(g::DiscretizedGrid{1}, coordinate::Float64)
    origcoord_to_grididx(g, (coordinate,))[1]
end
