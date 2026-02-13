struct InherentDiscreteGrid{D} <: Grid{D}
    Rs::NTuple{D,Int}
    origin::NTuple{D,Int}
    step::NTuple{D,Int}
    variablenames::NTuple{D,Symbol}
    base::NTuple{D,Int}
    indextable::Vector{Vector{Tuple{Symbol,Int}}}

    # Lookup table: lookup_table[variablename_index][bitnumber] -> (site_index, position_in_site)
    lookup_table::NTuple{D,Vector{Tuple{Int,Int}}}
    maxgrididx::NTuple{D,Int}
    site_radices::Vector{Vector{Int}}
    site_placevalues::Vector{Vector{Int}}
    sitedims::Vector{Int}

    function InherentDiscreteGrid{D}(
        Rs::NTuple{D,Int},
        origin::NTuple{D,Int},
        step::NTuple{D,Int},
        variablenames::NTuple{D,Symbol},
        base::NTuple{D,Int},
        indextable::Vector{Vector{Tuple{Symbol,Int}}}
    ) where {D}
        if !(D isa Int)
            throw(ArgumentError(lazy"Got dimension $D, which is not an Int."))
        end
        for d in 1:D
            if base[d] <= 1
                throw(ArgumentError(lazy"Got base[$d] = $(base[d]). base must be at least 2."))
            end
            if !(step[d] >= 1)
                throw(ArgumentError(lazy"Got step[$d] = $(step[d]). step must be at least 1."))
            end
            if !rangecheck_R(Rs[d]; base=base[d])
                throw(ArgumentError(lazy"Got Rs[$d] = $(Rs[d]) with base = $(base[d]). For all gridindices from 1 to base^R to fit into an Int, we must have base ^ R <= typemax(Int)"))
            end
        end
        if !allunique(variablenames)
            throw(ArgumentError(lazy"Got variablenames = $variablenames. variablenames must be unique."))
        end

        lookup_table = _build_lookup_table(Rs, indextable, variablenames)
        maxgrididx = map((b, R) -> b^R, base, Rs)
        site_radices = _build_site_radices(indextable, variablenames, base)
        site_placevalues = _build_site_placevalues(site_radices)
        sitedims = _build_site_dims(site_radices)

        return new{D}(Rs, origin, step, variablenames, base, indextable, lookup_table, maxgrididx, site_radices, site_placevalues, sitedims)
    end
end

# ============================================================================
# Helper/utility functions
# ============================================================================

_convert_to_scalar_if_possible(x) = x
_convert_to_scalar_if_possible(x::NTuple{1,T}) where {T} = first(x)

_to_tuple(::Val{d}, x::NTuple{d}) where {d} = x
_to_tuple(::Val{d}, x) where {d} = ntuple(i -> x, d)

_all_base_two(base::NTuple{D,Int}) where {D} = all(==(2), base)

function _normalize_indextable(indextable)
    normalized = Vector{Vector{Tuple{Symbol,Int}}}(undef, length(indextable))
    for (site_idx, site) in pairs(indextable)
        if !(site isa AbstractVector)
            throw(ArgumentError(lazy"Index table entry at site $site_idx must be an AbstractVector, got $(typeof(site))."))
        end

        normalized_site = Vector{Tuple{Symbol,Int}}(undef, length(site))
        for (position, qindex) in pairs(site)
            if !(qindex isa Tuple) || length(qindex) != 2
                throw(ArgumentError(lazy"Index table entry at site $site_idx position $position must be a tuple (Symbol, Int)."))
            end

            variablename, bitnumber = qindex
            if !(variablename isa Symbol)
                throw(ArgumentError(lazy"Index table entry at site $site_idx position $position has invalid variablename $(variablename). Expected Symbol."))
            end
            if !(bitnumber isa Integer)
                throw(ArgumentError(lazy"Index table entry at site $site_idx position $position has invalid bitnumber $(bitnumber). Expected Integer."))
            end

            normalized_site[position] = (variablename, Int(bitnumber))
        end

        normalized[site_idx] = normalized_site
    end

    return normalized
end

function _base_as_scalar_if_uniform(base::NTuple{D,Int}) where {D}
    if D == 0
        return base
    end
    return allequal(base) ? first(base) : base
end

default_step(::Val{D}) where D = ntuple(Returns(1), D)

default_origin(::Val{D}) where D = ntuple(Returns(1), D)

function _build_lookup_table(Rs::NTuple{D,Int}, indextable::Vector{Vector{Tuple{Symbol,Int}}}, variablenames::NTuple{D,Symbol}) where D
    lookup_table = ntuple(D) do d
        if Rs[d] < 0
            throw(ArgumentError(lazy"Got Rs[$d] = $(Rs[d]). Rs must be non-negative."))
        end
        Vector{Tuple{Int,Int}}(undef, Rs[d])
    end

    index_visited = [fill(false, Rs[d]) for d in 1:D]

    for (site_idx, quanticsindices) in pairs(indextable)
        for (pos_in_site, qindex) in pairs(quanticsindices)
            variablename, bitnumber = qindex
            var_idx = findfirst(==(variablename), variablenames)
            if isnothing(var_idx)
                throw(ArgumentError(lazy"Index table contains unknown index $qindex. Valid variablenames are $variablenames."))
            elseif bitnumber > Rs[var_idx]
                throw(ArgumentError(lazy"Index table contains quantics index $bitnumber of variable $variablename, but it must be smaller than or equal to the number of quantics indices for that variable, which is $(Rs[var_idx])."))
            elseif index_visited[var_idx][bitnumber]
                throw(ArgumentError(lazy"Index table contains quantics index $bitnumber of variable $variablename more than once."))
            end

            lookup_table[var_idx][bitnumber] = (site_idx, pos_in_site)
            index_visited[var_idx][bitnumber] = true
        end
    end

    for (var_idx, visited) in enumerate(index_visited)
        bitnumber = findfirst(==(false), visited)
        if !isnothing(bitnumber)
            throw(ArgumentError(lazy"Index table contains no site for quantics index $bitnumber of variable $(variablenames[var_idx])."))
        end
    end

    return lookup_table
end

function _build_site_radices(indextable::Vector{Vector{Tuple{Symbol,Int}}}, variablenames::NTuple{D,Symbol}, base::NTuple{D,Int}) where {D}
    var_index = Dict(zip(variablenames, 1:D))
    site_radices = Vector{Vector{Int}}(undef, length(indextable))
    for (i, site) in pairs(indextable)
        radices = Vector{Int}(undef, length(site))
        for (j, (variablename, _)) in pairs(site)
            radices[j] = base[var_index[variablename]]
        end
        site_radices[i] = radices
    end
    return site_radices
end

function _build_site_placevalues(site_radices::Vector{Vector{Int}})
    return map(site_radices) do radices
        len = length(radices)
        placevalues = Vector{Int}(undef, len)
        mult = 1
        for pos in len:-1:1
            placevalues[pos] = mult
            mult *= radices[pos]
        end
        placevalues
    end
end

function _build_site_dims(site_radices::Vector{Vector{Int}})
    return map(site_radices) do radices
        prod(radices; init=1)
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

function _build_indextable(variablenames::NTuple{D,Symbol}, Rs::NTuple{D,Int}, unfoldingscheme::Symbol) where D
    if !(unfoldingscheme in (:interleaved, :fused, :grouped))
        throw(ArgumentError(lazy"Got unfoldingscheme = $unfoldingscheme. Supported are :interleaved, :fused and :grouped."))
    end
    indextable = Vector{Tuple{Symbol,Int}}[]

    if unfoldingscheme === :grouped
        for d in 1:D
            for bitnumber in 1:Rs[d]
                qindex = (variablenames[d], bitnumber)
                push!(indextable, [qindex])
            end
        end
    else
        for bitnumber in 1:maximum(Rs; init=0)
            if unfoldingscheme === :interleaved
                _add_interleaved_indices!(indextable, variablenames, Rs, bitnumber)
            elseif unfoldingscheme === :fused
                _add_fused_indices!(indextable, variablenames, Rs, bitnumber)
            end
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

function _quantics_to_grididx_general(g::InherentDiscreteGrid{D}, quantics) where D
    return ntuple(D) do d
        grididx = 1
        R_d = g.Rs[d]

        for bitnumber in 1:R_d
            site_idx, pos_in_site = g.lookup_table[d][bitnumber]
            quantics_val = quantics[site_idx]
            zero_based = quantics_val - 1
            placevalue = g.site_placevalues[site_idx][pos_in_site]
            base_pos = g.site_radices[site_idx][pos_in_site]
            digit = (zero_based ÷ placevalue) % base_pos

            grididx += digit * g.base[d]^(R_d - bitnumber)
        end
        grididx
    end
end

function _quantics_to_grididx_base2(g::InherentDiscreteGrid{D}, quantics) where D
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

function _grididx_to_quantics_general!(result::Vector{Int}, g::InherentDiscreteGrid{D}, grididx::NTuple{D,Int}) where D
    @inbounds for d in 1:D
        zero_based_idx = grididx[d] - 1
        R_d = g.Rs[d]

        for bitnumber in 1:R_d
            site_idx, pos_in_site = g.lookup_table[d][bitnumber]
            base_d = g.base[d]
            bit_position = R_d - bitnumber
            digit = (zero_based_idx ÷ (base_d^bit_position)) % base_d

            placevalue = g.site_placevalues[site_idx][pos_in_site]
            result[site_idx] += digit * placevalue
        end
    end
end

function _grididx_to_quantics_base2!(result::Vector{Int}, g::InherentDiscreteGrid{D}, grididx::NTuple{D,Int}) where D
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

function InherentDiscreteGrid{D}(
    Rs,
    origin=default_origin(Val(D));
    unfoldingscheme=:fused,
    step=default_step(Val(D)),
    base=2
) where D
    Rs = _to_tuple(Val(D), Rs)
    origin = _to_tuple(Val(D), origin)
    step = _to_tuple(Val(D), step)
    base = _to_tuple(Val(D), base)
    variablenames = ntuple(Symbol, D)
    indextable = _build_indextable(variablenames, Rs, unfoldingscheme)
    InherentDiscreteGrid{D}(Rs, origin, step, variablenames, base, indextable)
end

function InherentDiscreteGrid(
    variablenames::NTuple{D,Symbol},
    indextable::AbstractVector;
    origin=default_origin(Val(D)),
    step=default_step(Val(D)),
    base=2
) where D
    indextable = _normalize_indextable(indextable)
    base = _to_tuple(Val(D), base)
    Rs = map(variablenames) do variablename
        count(index -> first(index) == variablename, Iterators.flatten(indextable))
    end

    return InherentDiscreteGrid{D}(Rs, origin, step, variablenames, base, indextable)
end

function InherentDiscreteGrid(
    variablenames::NTuple{D,Symbol},
    Rs::NTuple{D,Int};
    origin=default_origin(Val(D)),
    step=default_step(Val(D)),
    base=2,
    unfoldingscheme=:fused
) where {D}
    origin = _to_tuple(Val(D), origin)
    step = _to_tuple(Val(D), step)
    base = _to_tuple(Val(D), base)
    indextable = _build_indextable(variablenames, Rs, unfoldingscheme)
    return InherentDiscreteGrid{D}(Rs, origin, step, variablenames, base, indextable)
end

function InherentDiscreteGrid(Rs::NTuple{D,Int}; variablenames=ntuple(Symbol, D), kwargs...) where {D}
    return InherentDiscreteGrid(variablenames, Rs; kwargs...)
end

function InherentDiscreteGrid(R::Int, origin; kwargs...)
    return InherentDiscreteGrid{length(origin)}(R, origin; kwargs...)
end

# ============================================================================
# Basic property accessor functions
# ============================================================================

Base.ndims(::InherentDiscreteGrid{D}) where D = D

Base.length(g::InherentDiscreteGrid) = length(g.indextable)

grid_Rs(g::InherentDiscreteGrid) = g.Rs

grid_indextable(g::InherentDiscreteGrid) = g.indextable

grid_bases(g::InherentDiscreteGrid) = g.base

grid_base(g::InherentDiscreteGrid) = _base_as_scalar_if_uniform(g.base)

grid_variablenames(g::InherentDiscreteGrid) = g.variablenames

grid_step(g::InherentDiscreteGrid) = _convert_to_scalar_if_possible(g.step)

grid_origin(g::InherentDiscreteGrid) = _convert_to_scalar_if_possible(g.origin)

function sitedim(g::InherentDiscreteGrid, site::Int)::Int
    if !(site ∈ eachindex(g.indextable))
        throw(DomainError(site, lazy"Site index out of bounds [1, $(length(g.indextable))]."))
    end
    return g.sitedims[site]
end

# ============================================================================
# Grid coordinate functions
# ============================================================================

grid_min(g::InherentDiscreteGrid) = _convert_to_scalar_if_possible(g.origin)

grid_max(g::InherentDiscreteGrid) = _convert_to_scalar_if_possible(
    grid_origin(g) .+ grid_step(g) .* (g.base .^ g.Rs .- 1),
)

# ============================================================================
# Core conversion functions
# ============================================================================

function quantics_to_grididx(g::InherentDiscreteGrid{D}, quantics::AbstractVector{Int}) where D
    # TODO: add switch to turn off input validation
    if !(length(quantics) == length(g))
        throw(ArgumentError(lazy"Quantics vector must have length $(length(g.indextable)), got $(length(quantics))."))
    end

    for site in eachindex(quantics)
        if !(1 <= quantics[site] <= sitedim(g, site))
            throw(DomainError(quantics[site], lazy"Quantics value for site $site out of range 1:$(sitedim(g, site))."))
        end
    end

    result = if _all_base_two(g.base)
        _quantics_to_grididx_base2(g, quantics)
    else
        _quantics_to_grididx_general(g, quantics)
    end

    return _convert_to_scalar_if_possible(result)
end

function grididx_to_quantics(g::InherentDiscreteGrid{D}, grididx::Int) where D
    if D != 1
        throw(ArgumentError(lazy"grididx_to_quantics(grid, ::Integer) is only supported for 1D grids."))
    end
    grididx_tuple = _to_tuple(Val(D), grididx)
    return grididx_to_quantics(g, grididx_tuple)
end
function grididx_to_quantics(g::InherentDiscreteGrid{D}, grididx_tuple::NTuple{D,Int}) where D
    # TODO: add switch to turn off input validation
    for d in 1:D
        if !(1 <= grididx_tuple[d] <= g.maxgrididx[d])
            throw(DomainError(grididx_tuple[d], lazy"Grid index out of bounds [1, $(g.maxgrididx[d])]."))
        end
    end

    result = ones(Int, length(g.indextable))
    if _all_base_two(g.base)
        _grididx_to_quantics_base2!(result, g, grididx_tuple)
    else
        _grididx_to_quantics_general!(result, g, grididx_tuple)
    end
    return result
end

function grididx_to_origcoord(g::InherentDiscreteGrid{D}, grididx) where D
    grididx = _to_tuple(Val(D), grididx)
    for d in 1:D
        if !(grididx[d] ∈ 1:g.maxgrididx[d])
            throw(DomainError(grididx[d], lazy"Grid index out of bounds [1, $(g.maxgrididx[d])]."))
        end
    end

    res = grid_origin(g) .+ (grididx .- 1) .* grid_step(g)

    return _convert_to_scalar_if_possible(res)
end

function origcoord_to_grididx(g::InherentDiscreteGrid{D}, coordinate) where {D}
    coord_tuple = _to_tuple(Val(D), coordinate)
    bounds_lower = grid_min(g)
    bounds_upper = grid_max(g)
    # TODO: think about the correct bounds to use here
    for d in 1:D
        if !(bounds_lower[d] <= coord_tuple[d] <= bounds_upper[d])
            throw(DomainError(coord_tuple[d], lazy"Coordinate out of bounds [$(bounds_lower[d]), $(bounds_upper[d])]."))
        end
    end

    discrete_idx = div.(coord_tuple .- grid_min(g), grid_step(g)) .+ 1
    discrete_idx = clamp.(discrete_idx, 1, g.base .^ g.Rs)

    return _convert_to_scalar_if_possible(discrete_idx)
end

function origcoord_to_quantics(g::InherentDiscreteGrid, coordinate)
    grididx_to_quantics(g, origcoord_to_grididx(g, coordinate))
end

function quantics_to_origcoord(g::InherentDiscreteGrid, quantics)
    grididx_to_origcoord(g, quantics_to_grididx(g, quantics))
end

# ============================================================================
# Other utility functions
# ============================================================================

function localdimensions(g::InherentDiscreteGrid)::Vector{Int}
    return copy(g.sitedims)
end
