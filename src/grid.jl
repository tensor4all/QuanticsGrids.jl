
# Grid for d-dimensional space
abstract type Grid{d} end


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
    return (index .- 1) .* grid_step(g) .+ grid_min(g)
end

"""
grid index => quantics
"""
function grididx_to_quantics(g::Grid{d}, grididx::NTuple{d,Int}) where {d}
    if g.unfoldingscheme == UnfoldingSchemes.fused
        return index_to_quantics_fused(grididx, numdigits=g.R, base=g.base)
    else
        return fused_to_interleaved(
            index_to_quantics_fused(grididx, numdigits=g.R, base=g.base),
            d,
            base=g.base
        )
    end
end

"""
original coordinate => quantics
"""
function origcoord_to_quantics(g::Grid{d}, coordinate::NTuple{d,T}) where {d,T}
    return grididx_to_quantics(g, origcoord_to_grididx(g, coordinate))
end


"""
quantics => grid index
"""
function quantics_to_grididx(g::Grid{d}, bitlist) where {d}
    return quantics_to_index(
        bitlist;
        base=g.base,
        dims=Val(d),
        unfoldingscheme=g.unfoldingscheme,
    )
end

"""
quantics => original coordinate system
"""
function quantics_to_origcoord(g::Grid{d}, bitlist) where {d}
    return grididx_to_origcoord(g, quantics_to_grididx(g, bitlist))
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
    unfoldingscheme::UnfoldingSchemes.UnfoldingScheme
    step::NTuple{d,Int}

    function InherentDiscreteGrid{d}(
        R::Int,
        origin::NTuple{d,Int};
        base::Integer=2,
        unfoldingscheme::UnfoldingSchemes.UnfoldingScheme=UnfoldingSchemes.fused,
        step::NTuple{d,Int}=ntuple(i -> 1, d),
    ) where {d}
        new(R, origin, base, unfoldingscheme, step)
    end
end

grid_min(grid::InherentDiscreteGrid) = grid.origin
grid_step(grid::InherentDiscreteGrid{d}) where {d} = grid.step

"""
Create a grid for inherently discrete data with origin at 1
"""
function InherentDiscreteGrid{d}(
    R::Int;
    base::Integer=2,
    step::NTuple{d,Int}=ntuple(i -> 1, d),
    unfoldingscheme::UnfoldingSchemes.UnfoldingScheme=UnfoldingSchemes.fused
) where {d}
    InherentDiscreteGrid{d}(
        R, ntuple(i -> 1, d); base=base, unfoldingscheme=unfoldingscheme, step=step
    )
end


"""
Create a grid for inherently discrete data with origin at an arbitrary point
"""
function InherentDiscreteGrid{d}(
    R::Int, origin::Int;
    base=2, unfoldingscheme::UnfoldingSchemes.UnfoldingScheme=UnfoldingSchemes.fused,
    step::NTuple{d,Int}=ntuple(i -> 1, d),
) where {d}
    InherentDiscreteGrid{d}(
        R, ntuple(i -> origin, d); base=base, unfoldingscheme=unfoldingscheme, step=step
    )
end


"""
Convert a coordinate in the original coordinate system to the corresponding grid index
"""
function origcoord_to_grididx(
    g::InherentDiscreteGrid,
    coordinate::Union{Int,NTuple{N,Int}},
) where {N}
    return div.(coordinate .- grid_min(g), grid_step(g)) .+ 1
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
    grid_min::NTuple{d,Float64}
    grid_max::NTuple{d,Float64}
    base::Int
    unfoldingscheme::UnfoldingSchemes.UnfoldingScheme
    includeendpoint::Bool

    function DiscretizedGrid{d}(
        R::Int, grid_min, grid_max;
        base::Integer=2,
        unfoldingscheme::UnfoldingSchemes.UnfoldingScheme=UnfoldingSchemes.fused,
        includeendpoint::Bool=false
    ) where {d}
        return new(R, grid_min, grid_max, base, unfoldingscheme, includeendpoint)
    end
end

grid_min(g::DiscretizedGrid) = g.grid_min
grid_max(g::DiscretizedGrid) = g.grid_max
grid_step(g::DiscretizedGrid{d}) where {d} = g.includeendpoint ? (g.grid_max .- g.grid_min) ./ (g.base^g.R - 1) : (g.grid_max .- g.grid_min) ./ (g.base^g.R)


function DiscretizedGrid{d}(
    R::Int;
    base=2,
    unfoldingscheme::UnfoldingSchemes.UnfoldingScheme=UnfoldingSchemes.fused
) where {d}
    return DiscretizedGrid{d}(
        R, ntuple(i -> 0.0, d), ntuple(i -> 1.0, d);
        base=base, unfoldingscheme=unfoldingscheme
    )
end


"""
Convert a coordinate in the original coordinate system to the corresponding grid index
"""
function origcoord_to_grididx(g::DiscretizedGrid, coordinate::NTuple{N,Float64}) where {N}
    if g.includeendpoint
        all(grid_min(g) .<= coordinate .<= grid_max(g)) ||
            error("Bound Error: $(coordinate), min=$(grid_min(g)), max=$(grid_max(g))")
    else
        all(grid_min(g) .<= coordinate .< grid_max(g)) ||
            error("Bound Error: $(coordinate), min=$(grid_min(g)), max=$(grid_max(g))")
    end
    return ((coordinate .- grid_min(g)) ./ grid_step(g) .+ 1) .|> floor .|> Int
end

function origcoord_to_grididx(g::DiscretizedGrid{1}, coordinate::Float64)
    origcoord_to_grididx(g, (coordinate,))[1]
end
