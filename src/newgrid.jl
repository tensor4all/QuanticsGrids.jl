@kwdef struct QuanticsIndex
    variablename::Symbol
    bitnumber::Int

    function QuanticsIndex(variablename::Symbol, bitnumber::Int)
        @assert bitnumber > 0
        return new(variablename, bitnumber)
    end
end

@kwdef struct NewDiscretizedGrid{D}
    Rs::NTuple{D,Int}
    lower_bound::NTuple{D,Float64}
    upper_bound::NTuple{D,Float64}
    variablenames::NTuple{D,Symbol}
    base::Int
    indextable::Vector{Vector{QuanticsIndex}}
    # Lookup table: (variablename_index, bitnumber) -> (site_index, position_in_site)
    lookup_table::Dict{Tuple{Int,Int},Tuple{Int,Int}}

    function NewDiscretizedGrid{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable) where {D}
        @assert all(>=(0), Rs)
        @assert all(lower_bound .< upper_bound)
        @assert base > 1
        @assert all(R -> rangecheck_R(R; base), Rs)

        # Build lookup table for fast access
        lookup_table = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
        for (site_idx, quanticsindices) in pairs(indextable)
            for (pos_in_site, qindex) in pairs(quanticsindices)
                var_idx = findfirst(==(qindex.variablename), variablenames)
                if !isnothing(var_idx)
                    lookup_table[var_idx, qindex.bitnumber] = (site_idx, pos_in_site)
                end
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

function NewDiscretizedGrid(Rs::NTuple{D,Int}; lower_bound=default_lower_bound(Val(D)), upper_bound=default_upper_bound(Val(D)), base=2) where {D}
    variablenames = ntuple(Symbol, D)
    indextable = Vector{QuanticsIndex}[]

    for bitnumber in 1:maximum(Rs)
        for d in 1:D
            bitnumber ∈ 1:Rs[d] || continue
            variablename = Symbol(d)
            qindex = QuanticsIndex(; variablename, bitnumber)
            push!(indextable, [qindex])
        end
    end

    return NewDiscretizedGrid{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable)
end

function NewDiscretizedGrid{D}(R::Int; lower_bound=default_lower_bound(Val(D)), upper_bound=default_upper_bound(Val(D)), base=2) where {D}
    return NewDiscretizedGrid(ntuple(Returns(R), D); lower_bound, upper_bound, base)
end

function NewDiscretizedGrid(variablenames::NTuple{D,Symbol}, indextable::Vector{Vector{QuanticsIndex}}; lower_bound=default_lower_bound(Val(D)), upper_bound=default_upper_bound(Val(D)), base=2) where D
    @assert all(Iterators.flatten(indextable)) do index
        index.variablename ∈ variablenames
    end

    Rs = Tuple(map(variablenames) do variablename
        count(index -> index.variablename == variablename, Iterators.flatten(indextable))
    end)
    NewDiscretizedGrid{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable)
end

function NewDiscretizedGrid(variablenames::NTuple{D,Symbol}, tupletable::Vector{Vector{Tuple{Symbol,Int}}}; lower_bound=default_lower_bound(Val(D)), upper_bound=default_upper_bound(Val(D)), base=2) where D
    indextable = map(tupletable) do tuplevec
        map(tuplevec) do tuple
            variablename, bitnumber = tuple
            QuanticsIndex(; variablename, bitnumber)
        end
    end
    NewDiscretizedGrid(variablenames, indextable; lower_bound, upper_bound, base)
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
    ntuple(D) do d
        grididx::Int = 1
        R_d = g.Rs[d]

        # Process each bitnumber for this dimension
        for bitnumber in 1:R_d
            # Look up site and position using precomputed lookup table
            site_idx, pos_in_site = g.lookup_table[d, bitnumber]
            quantics_val = quantics[site_idx]
            site_len = length(g.indextable[site_idx])

            # Extract the digit at this position without allocation
            # Convert 1-based quantics value to 0-based, then extract digit
            temp = quantics_val - 1
            for _ in 1:(site_len-pos_in_site)
                temp = div(temp, base)
            end
            digit = temp % base

            # Add contribution to grid index
            grididx += digit * base^(R_d - bitnumber)
            # end
        end
        grididx
    end
end

function grididx_to_quantics(g::NewDiscretizedGrid{D}, grididx) where D
    @assert length(grididx) == D
    digs = [reverse(digits(grididx[d] - 1; base=g.base, pad=g.Rs[d])) for d in 1:D]

    map(g.indextable) do quanticsindices
        site_value = 0
        for (idx, quanticsindex) in enumerate(quanticsindices)
            d = findfirst(==(quanticsindex.variablename), g.variablenames)
            site_value += digs[d][quanticsindex.bitnumber] * g.base^(length(quanticsindices) - idx)
        end
        site_value + 1
    end
end
