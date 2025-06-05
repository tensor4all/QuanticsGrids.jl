@kwdef struct QuanticsIndex
    variablename::Symbol
    bitnumber::Int

    function QuanticsIndex(variablename::Symbol, bitnumber::Int)
        @assert bitnumber > 0
        return new(variablename, bitnumber)
    end
end

struct NewDiscretizedGrid{D}
    Rs::NTuple{D,Int}
    lower_bound::NTuple{D,Float64}
    upper_bound::NTuple{D,Float64}
    variablenames::NTuple{D,Symbol}
    base::Int
    indextable::Vector{Vector{QuanticsIndex}}
    # Lookup table: lookup_table[variablename_index][bitnumber] -> (site_index, position_in_site)
    lookup_table::NTuple{D,Vector{Tuple{Int,Int}}}

    function NewDiscretizedGrid{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable, includeendpoint) where {D}
        @assert all(>=(0), Rs)
        @assert all(lower_bound .< upper_bound)
        @assert base > 1
        @assert all(R -> rangecheck_R(R; base), Rs)

        # Build lookup table for fast access
        lookup_table = ntuple(d -> Vector{Tuple{Int,Int}}(undef, Rs[d]), D)

        for (site_idx, quanticsindices) in pairs(indextable)
            for (pos_in_site, qindex) in pairs(quanticsindices)
                var_idx = findfirst(==(qindex.variablename), variablenames)
                if !isnothing(var_idx)
                    lookup_table[var_idx][qindex.bitnumber] = (site_idx, pos_in_site)
                end
            end
        end

        if includeendpoint isa Bool
            includeendpoint = ntuple(Returns(includeendpoint), D)
        end

        if lower_bound isa Number
            lower_bound = ntuple(Returns(lower_bound), D)
        end
        if upper_bound isa Number
            upper_bound = ntuple(Returns(upper_bound), D)
        end

        upper_bound = ntuple(D) do d
            if includeendpoint[d]
                upper_bound[d] + (upper_bound[d] - lower_bound[d]) / (base^Rs[d] - 1)
            else
                upper_bound[d]
            end
        end

        new{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable, lookup_table)
    end
end

function rangecheck_R(R; base=2)
    # For all gridindices from 1 to base^R to fit into an Int64,
    # we must have base ^ R <= typemax(Int)
    R * log(base) <= log(typemax(Int))
end

function NewDiscretizedGrid(Rs::NTuple{D,Int}; lower_bound=default_lower_bound(Val(D)), upper_bound=default_upper_bound(Val(D)), base=2, unfoldingscheme=:interleaved, includeendpoint=false) where {D}
    variablenames = ntuple(Symbol, D)
    indextable = Vector{QuanticsIndex}[]

    for bitnumber in 1:maximum(Rs)
        if unfoldingscheme == :interleaved
            for d in 1:D
                bitnumber ∈ 1:Rs[d] || continue
                variablename = Symbol(d)
                qindex = QuanticsIndex(; variablename, bitnumber)
                push!(indextable, [qindex])
            end
        elseif unfoldingscheme == :fused
            indices_bitnumber = Vector{QuanticsIndex}()
            for d in 1:D
                bitnumber ∈ 1:Rs[d] || continue
                variablename = Symbol(d)
                qindex = QuanticsIndex(; variablename, bitnumber)
                push!(indices_bitnumber, qindex)
            end
            push!(indextable, indices_bitnumber)
        else
            error(lazy"""Unfolding scheme $unfoldingscheme not supported. Use :interleaved or :fused.
            If you need a different scheme, please use the NewDiscretizedGrid(variablenames::NTuple{D,Symbol}, indextable::Vector{Vector{QuanticsIndex}}) constructor.""")
        end
    end

    return NewDiscretizedGrid{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable, includeendpoint)
end

function NewDiscretizedGrid{D}(R::Int, lower_bound=default_lower_bound(Val(D)), upper_bound=default_upper_bound(Val(D)); kwargs...) where {D}
    return NewDiscretizedGrid(ntuple(Returns(R), D); lower_bound, upper_bound, kwargs...)
end

function NewDiscretizedGrid(variablenames::NTuple{D,Symbol}, indextable::Vector{Vector{QuanticsIndex}}; lower_bound=default_lower_bound(Val(D)), upper_bound=default_upper_bound(Val(D)), base=2, includeendpoint=false) where D
    @assert all(Iterators.flatten(indextable)) do index
        index.variablename ∈ variablenames
    end

    Rs = Tuple(map(variablenames) do variablename
        count(index -> index.variablename == variablename, Iterators.flatten(indextable))
    end)
    NewDiscretizedGrid{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable, includeendpoint)
end

function NewDiscretizedGrid(variablenames::NTuple{D,Symbol}, tupletable::Vector{Vector{Tuple{Symbol,Int}}}; kwargs...) where D
    indextable = map(tupletable) do tuplevec
        map(tuplevec) do tuple
            variablename, bitnumber = tuple
            QuanticsIndex(; variablename, bitnumber)
        end
    end
    NewDiscretizedGrid(variablenames, indextable; kwargs...)
end

default_lower_bound(::Val{D}) where D = ntuple(Returns(0.0), D)
default_upper_bound(::Val{D}) where D = ntuple(Returns(1.0), D)

Base.ndims(::NewDiscretizedGrid{D}) where D = D
Base.length(g::NewDiscretizedGrid) = length(g.indextable)

function sitedim(g::NewDiscretizedGrid, site)
    @assert site ∈ eachindex(g.indextable)
    return g.base^length(g.indextable[site])
end

function quantics_to_grididx(g::NewDiscretizedGrid{D}, quantics) where D
    @assert length(quantics) == length(g)
    @assert all(site -> quantics[site] ∈ 1:sitedim(g, site), eachindex(quantics))

    base = g.base

    result = if base == 2
        _quantics_to_grididx_base2(g, quantics)
    else
        ntuple(D) do d
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
    _convert_to_scalar_if_possible(result)
end

function _quantics_to_grididx_base2(g::NewDiscretizedGrid{D}, quantics) where D
    ntuple(D) do d
        grididx = 0
        R_d = g.Rs[d]

        for bitnumber in 1:R_d
            site_idx, pos_in_site = g.lookup_table[d][bitnumber]
            bit_position = length(g.indextable[site_idx]) - pos_in_site
            digit = (quantics[site_idx] - 1) >> bit_position & 1
            grididx |= digit << (R_d - bitnumber)
        end
        grididx + 1  # Convert back to 1-based indexing
    end
end

"""
    grididx_to_quantics(g::NewDiscretizedGrid{D}, grididx::NTuple{D,Int}) where D

Convert grid indices to quantics representation.
"""
function grididx_to_quantics(g::NewDiscretizedGrid{D}, grididx::NTuple{D,Int}) where D
    @assert length(grididx) == D "Grid index must have dimension $D, got $(length(grididx))"
    @assert all(1 ≤ grididx[d] ≤ g.base^g.Rs[d] for d in 1:D) "Grid index out of bounds"

    base = g.base

    result = ones(Int, length(g.indextable))

    if base == 2
        _grididx_to_quantics_base2!(result, g, grididx)
    else
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

    return result
end

grididx_to_quantics(g::NewDiscretizedGrid{1}, grididx::Int) = grididx_to_quantics(g, (grididx,))

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

"""
    grid_min(g::NewDiscretizedGrid)

Returns the grid point with minimal coordinate values.
This is equivalent to the lower bound.
The return value is scalar for 1D grids, and a `Tuple` otherwise.
"""
grid_min(g::NewDiscretizedGrid) = _convert_to_scalar_if_possible(g.lower_bound)

"""
    grid_max(g::NewDiscretizedGrid)

Returns the grid point with maximal coordinate values.
Since NewDiscretizedGrid doesn't have includeendpoint, this is always upper_bound - grid_step.
The return value is scalar for 1D grids, and a `Tuple` otherwise.
"""
grid_max(g::NewDiscretizedGrid) = _convert_to_scalar_if_possible(g.upper_bound .- grid_step(g))

"""
    grid_step(g::NewDiscretizedGrid)

Returns the distance between adjacent grid points in each dimension.
The return value is scalar for 1D grids, and a `Tuple` otherwise.
"""
grid_step(g::NewDiscretizedGrid) = _convert_to_scalar_if_possible(
    (upper_bound(g) .- lower_bound(g)) ./ (g.base .^ g.Rs),
)

"""
    upper_bound(g::NewDiscretizedGrid)

Returns the upper bound of the grid coordinates.
The return value is scalar for 1D grids, and a `Tuple` otherwise.
"""
upper_bound(g::NewDiscretizedGrid) = _convert_to_scalar_if_possible(g.upper_bound)

"""
    lower_bound(g::NewDiscretizedGrid)

Returns the lower bound of the grid coordinates.
The return value is scalar for 1D grids, and a `Tuple` otherwise.
"""
lower_bound(g::NewDiscretizedGrid) = _convert_to_scalar_if_possible(g.lower_bound)

"""
Convert a grid index to the corresponding coordinate in the original coordinate system
"""
function grididx_to_origcoord(g::NewDiscretizedGrid{D}, index::NTuple{D,Int}) where {D}
    all(1 .<= index .<= (g.base .^ g.Rs)) || error(lazy"Grid-index $index out of bounds [1, $(g.base .^ g.Rs)]")
    return _convert_to_scalar_if_possible((index .- 1) .* grid_step(g) .+ grid_min(g))
end

"""
Convert a grid index to the corresponding coordinate in the original coordinate system (1D version)
"""
grididx_to_origcoord(g::NewDiscretizedGrid{1}, index::Int) = grididx_to_origcoord(g, (index,))

"""
Convert a coordinate in the original coordinate system to the corresponding grid index
"""
function origcoord_to_grididx(g::NewDiscretizedGrid{D}, coordinate::NTuple{D,Float64}) where {D}
    if !all(lower_bound(g) .<= coordinate .<= upper_bound(g))
        throw(BoundsError(g, coordinate))
    end
    # Calculate raw index and clamp to valid range
    min_ = grid_min(g)
    step_ = grid_step(g)
    raw_indices = @. round(Int, (coordinate - min_) / step_ + 1)
    max_indices = g.base .^ g.Rs
    min_indices = ntuple(Returns(1), D)
    clamped_indices = clamp.(raw_indices, min_indices, max_indices)
    return _convert_to_scalar_if_possible(clamped_indices)
end

"""
Convert a coordinate in the original coordinate system to the corresponding grid index (1D version)
"""
function origcoord_to_grididx(g::NewDiscretizedGrid{1}, coordinate::Float64)
    origcoord_to_grididx(g, (coordinate,))
end

"""
Convert original coordinate to quantics representation
"""
function origcoord_to_quantics(g::NewDiscretizedGrid{d}, coordinate) where {d}
    return grididx_to_quantics(g, origcoord_to_grididx(g, _to_tuple(Val(d), coordinate)))
end

"""
Convert quantics representation to original coordinate
"""
function quantics_to_origcoord(g::NewDiscretizedGrid{d}, quantics) where {d}
    return _convert_to_scalar_if_possible(
        grididx_to_origcoord(g, quantics_to_grididx(g, quantics)),
    )
end

function localdimensions(g::NewDiscretizedGrid)
    return g.base .^ length.(g.indextable)
end

"""
Make a wrapper function that takes a bitlist as input
"""
function quanticsfunction(::Type{T}, g::NewDiscretizedGrid, f::F)::Function where {T,F<:Function}
    function _f(bitlist)::T
        return f(quantics_to_origcoord(g, bitlist)...)
    end
    return _f
end

grid_origin(g::NewDiscretizedGrid) = _convert_to_scalar_if_possible(
    lower_bound(g)
)

