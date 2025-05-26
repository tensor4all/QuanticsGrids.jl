using ReTestItems

@kwdef struct QuanticsIndex
    variablename::Symbol
    bitnumber::Int
end

@kwdef struct NewDiscretizedGrid{D}
    Rs::NTuple{D,Int}
    lower_bound::NTuple{D,Float64}
    upper_bound::NTuple{D,Float64}
    variablenames::NTuple{D,Symbol}
    base::Int
    indextable::Vector{Vector{QuanticsIndex}}

    function NewDiscretizedGrid{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable) where {D}
        @assert all(>=(0), Rs)
        @assert all(lower_bound .< upper_bound)
        @assert base > 1
        @assert all(R -> rangecheck_R(R; base), Rs)

        new{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable)
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

    return NewDiscretizedGrid{D}(; Rs, lower_bound, upper_bound, variablenames, base, indextable)
end

function NewDiscretizedGrid{D}(R::Int; lower_bound=default_lower_bound(Val(D)), upper_bound=default_upper_bound(Val(D)), base=2) where {D}
    return NewDiscretizedGrid(ntuple(Returns(R), D); lower_bound, upper_bound, base)
end

default_lower_bound(::Val{D}) where D = ntuple(Returns(0.0), D)
default_upper_bound(::Val{D}) where D = ntuple(Returns(1.0), D)

Base.ndims(::NewDiscretizedGrid{D}) where D = D
Base.length(g::NewDiscretizedGrid) = sum(g.Rs)

function sitedim(g::NewDiscretizedGrid, site)
    @assert site ∈ eachindex(g.indextable)
    return g.base^length(g.indextable[site])
end

function quantics_to_grididx(g::NewDiscretizedGrid{D}, quantics) where D
    @assert length(quantics) == length(g)
    @assert all(site -> quantics[site] ∈ 1:sitedim(g, site), eachindex(quantics))

    ntuple(D) do d
        variablename = g.variablenames[d]
        grididx = 1
        for (site, quanticsindices) in pairs(g.indextable)
            idx = findfirst(quanticsindex -> quanticsindex.variablename == variablename, quanticsindices)
            isnothing(idx) && continue
            bitnumber = quanticsindices[idx].bitnumber
            grididx += (quantics[site] - 1) * g.base^(g.Rs[d] - bitnumber)
        end
        grididx
    end
end

function grididx_to_quantics(g::NewDiscretizedGrid{D}, grididx) where D
    @assert length(grididx) == D
    digits_d = [reverse(digits(grididx[d] - 1; base=g.base, pad=g.Rs[d])) for d in 1:D]

    map(g.indextable) do quanticsindices
        for quanticsindex in quanticsindices
            d = findfirst(==(quanticsindex.variablename), g.variablenames)
            return digits_d[d][quanticsindex.bitnumber] + 1
        end
    end
end
