export NewInherentDiscreteGrid

struct NewInherentDiscreteGrid{D}
    Rs::NTuple{D,Int}
    origin::NTuple{D,Int}
    step::NTuple{D,Int}
    variablenames::NTuple{D,Symbol}
    base::Int
    indextable::Vector{Vector{Tuple{Symbol,Int}}}

    # Lookup table: lookup_table[variablename_index][bitnumber] -> (site_index, position_in_site)
    lookup_table::NTuple{D,Vector{Tuple{Int,Int}}}

    function NewInherentDiscreteGrid{D}(
        Rs::NTuple{D,Int},
        origin::NTuple{D,Int},
        step::NTuple{D,Int},
        variablenames::NTuple{D,Symbol},
        base::Int,
        indextable::Vector{Vector{Tuple{Symbol,Int}}}
    ) where {D}
        @assert all(>=(0), Rs)
        @assert all(>=(1), step)
        @assert base > 1
        @assert all(R -> rangecheck_R(R; base), Rs)
        @assert allunique(variablenames)

        lookup_table = _build_lookup_table(Rs, indextable, variablenames, Val(D))

        # origin = _to_tuple(Val(D), origin)
        # step = _to_tuple(Val(D), step)

        return new{D}(Rs, origin, step, variablenames, base, indextable, lookup_table)
    end
end

# ============================================================================
# Helper/utility functions
# ============================================================================

default_step(::Val{D}) where D = ntuple(Returns(1), D)

default_origin(::Val{D}) where D = ntuple(Returns(1), D)

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

function _build_indextable(variablenames::NTuple{D,Symbol}, Rs::NTuple{D,Int}, unfoldingscheme::Symbol) where D
    @assert unfoldingscheme in (:interleaved, :fused)
    indextable = Vector{Tuple{Symbol,Int}}[]

    for bitnumber in 1:maximum(Rs)
        if unfoldingscheme === :interleaved
            _add_interleaved_indices!(indextable, variablenames, Rs, bitnumber)
        elseif unfoldingscheme === :fused
            _add_fused_indices!(indextable, variablenames, Rs, bitnumber)
        end
    end

    return indextable
end

function _add_interleaved_indices!(indextable, variablenames::NTuple{D,Symbol}, Rs::NTuple{D,Int}, bitnumber) where D
    for d in 1:D
        bitnumber ∈ 1:Rs[d] || continue
        qindex = (variablenames[d], bitnumber)
        push!(indextable, [qindex])
    end
end

function _add_fused_indices!(indextable, variablenames::NTuple{D,Symbol}, Rs::NTuple{D,Int}, bitnumber) where D
    indices_bitnumber = Tuple{Symbol,Int}[]
    # Add dimensions in reverse order to match DiscretizedGrid convention
    # where the first dimension varies fastest in fused quantics
    for d in D:-1:1
        bitnumber ∈ 1:Rs[d] || continue
        qindex = (variablenames[d], bitnumber)
        push!(indices_bitnumber, qindex)
    end
    if !isempty(indices_bitnumber)
        push!(indextable, indices_bitnumber)
    end
end

function _quantics_to_grididx_general(g::NewInherentDiscreteGrid{D}, quantics) where D
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

function _quantics_to_grididx_base2(g::NewInherentDiscreteGrid{D}, quantics) where D
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

function _grididx_to_quantics_general!(result::Vector{Int}, g::NewInherentDiscreteGrid{D}, grididx::NTuple{D,Int}) where D
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

function _grididx_to_quantics_base2!(result::Vector{Int}, g::NewInherentDiscreteGrid{D}, grididx::NTuple{D,Int}) where D
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

# ============================================================================
# Constructors
# ============================================================================

function NewInherentDiscreteGrid{D}(
    Rs,
    origin=default_origin(Val(D));
    unfoldingscheme=:fused,
    step=default_step(Val(D)),
    base=2
) where D
    Rs = _to_tuple(Val(D), Rs)
    origin = _to_tuple(Val(D), origin)
    step = _to_tuple(Val(D), step)
    variablenames = ntuple(Symbol, D)
    indextable = _build_indextable(variablenames, Rs, unfoldingscheme)
    NewInherentDiscreteGrid{D}(Rs, origin, step, variablenames, base, indextable)
end

function NewInherentDiscreteGrid(
    variablenames::NTuple{D,Symbol},
    indextable::Vector{Vector{Tuple{Symbol,Int}}};
    origin=default_origin(Val(D)),
    step=default_step(Val(D)),
    base=2
) where D
    @assert all(Iterators.flatten(indextable)) do index
        first(index) ∈ variablenames
    end

    Rs = Tuple(map(variablenames) do variablename
        count(index -> first(index) == variablename, Iterators.flatten(indextable))
    end)

    return NewInherentDiscreteGrid{D}(Rs, origin, step, variablenames, base, indextable)
end

function NewInherentDiscreteGrid(
    variablenames::NTuple{D,Symbol},
    Rs::NTuple{D,Int};
    origin=default_origin(Val(D)),
    step=default_step(Val(D)),
    base=2,
    unfoldingscheme=:fused
) where {D}
    origin = _to_tuple(Val(D), origin)
    step = _to_tuple(Val(D), step)
    indextable = _build_indextable(variablenames, Rs, unfoldingscheme)
    return NewInherentDiscreteGrid{D}(Rs, origin, step, variablenames, base, indextable)
end

function NewInherentDiscreteGrid(Rs::NTuple{D,Int}; kwargs...) where {D}
    variablenames = ntuple(Symbol, D)
    return NewInherentDiscreteGrid(variablenames, Rs; kwargs...)
end

function NewInherentDiscreteGrid(R::Int, origin; kwargs...)
    return NewInherentDiscreteGrid{1}(R, origin; kwargs...)
end

# ============================================================================
# Basic property accessor functions
# ============================================================================

Base.ndims(::NewInherentDiscreteGrid{D}) where D = D

Base.length(g::NewInherentDiscreteGrid) = length(g.indextable)

grid_Rs(g::NewInherentDiscreteGrid) = g.Rs

grid_indextable(g::NewInherentDiscreteGrid) = g.indextable

grid_base(g::NewInherentDiscreteGrid) = g.base

grid_variablenames(g::NewInherentDiscreteGrid) = g.variablenames

grid_step(g::NewInherentDiscreteGrid) = _convert_to_scalar_if_possible(g.step)

grid_origin(g::NewInherentDiscreteGrid) = _convert_to_scalar_if_possible(g.origin)

function sitedim(g::NewInherentDiscreteGrid, site::Int)::Int
    @assert site ∈ eachindex(g.indextable)
    return g.base^length(g.indextable[site])
end

# ============================================================================
# Grid coordinate functions
# ============================================================================

grid_min(g::NewInherentDiscreteGrid) = _convert_to_scalar_if_possible(g.origin)

grid_max(g::NewInherentDiscreteGrid) = _convert_to_scalar_if_possible(
    grid_origin(g) .+ grid_step(g) .* (g.base .^ g.Rs .- 1),
)

# ============================================================================
# Core conversion functions
# ============================================================================

function quantics_to_grididx(g::NewInherentDiscreteGrid{D}, quantics::AbstractVector{Int}) where D
    @assert length(quantics) == length(g)
    @assert all(site -> quantics[site] ∈ 1:sitedim(g, site), eachindex(quantics))

    result = if g.base == 2
        _quantics_to_grididx_base2(g, quantics)
    else
        _quantics_to_grididx_general(g, quantics)
    end

    return _convert_to_scalar_if_possible(result)
end

function grididx_to_quantics(g::NewInherentDiscreteGrid{D}, grididx) where D
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

function grididx_to_origcoord(g::NewInherentDiscreteGrid{D}, grididx) where D
    grididx = _to_tuple(Val(D), grididx)
    @assert all(1 .<= grididx .<= (g.base .^ g.Rs)) lazy"Grid-index $grididx out of bounds [1, $(g.base .^ g.Rs)]"

    res = grid_origin(g) .+ (grididx .- 1) .* grid_step(g)

    return _convert_to_scalar_if_possible(res)
end

function origcoord_to_grididx(g::NewInherentDiscreteGrid{D}, coordinate) where {D}
    coord_tuple = _to_tuple(Val(D), coordinate)
    bounds_lower = grid_min(g)
    bounds_upper = grid_max(g)
    # TODO: think about the correct bounds to use here
    @assert all(bounds_lower .<= coord_tuple .<= bounds_upper) "Coordinate $coord_tuple out of bounds [$bounds_lower, $bounds_upper]"

    discrete_idx = div.(coord_tuple .- grid_min(g), grid_step(g)) .+ 1
    discrete_idx = clamp.(discrete_idx, 1, g.base .^ g.Rs)

    return _convert_to_scalar_if_possible(discrete_idx)
end

function origcoord_to_quantics(g::NewInherentDiscreteGrid, coordinate)
    grididx_to_quantics(g, origcoord_to_grididx(g, coordinate))
end

function quantics_to_origcoord(g::NewInherentDiscreteGrid, quantics)
    grididx_to_origcoord(g, quantics_to_grididx(g, quantics))
end

# ============================================================================
# Other utility functions
# ============================================================================

function localdimensions(g::NewInherentDiscreteGrid)::Vector{Int}
    return g.base .^ length.(g.indextable)
end
