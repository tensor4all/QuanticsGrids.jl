function _checkdigits(base::Integer, digitlist)
    maximum(digitlist) <= base || error("maximum(digitlist) <= base")
    minimum(digitlist) >= 0 || error("minimum(digitlist) >= 0")
end

"""
    fuse_dimensions(digitlists...)

Fuse d digitlists that represent a quantics index into a digitlist where each bit
has dimension base^d. This fuses legs for different dimensions that have equal length
scale (see QTCI paper).

Inverse of [`unfuse_dimensions`](@ref).
"""
function fuse_dimensions(digitlists...; base::Integer = 2)
    _checkdigits.(base, digitlists)
    result = ones(Int, length(digitlists[1]))
    return fuse_dimensions!(result, digitlists...; base = base)
end



function fuse_dimensions!(fused::AbstractArray{<:Integer}, digitlists...; base::Integer = 2)
    p = 1
    fused .= 1
    for d in eachindex(digitlists)
        @. fused += (digitlists[d] - 1) * p
        p *= base
    end
    return fused
end


"""
    unfuse_dimensions([base=Val(2)], digitlist, d)

Unfuse up a fused digitlist with bits of dimension base^d into d digitlists where each bit has dimension `base`.
Inverse of [`fuse_dimensions`](@ref).
"""
function unfuse_dimensions(digitlist, d; base::Integer = 2)
    result = [zeros(Int, length(digitlist)) for _ = 1:d]
    return unfuse_dimensions!(result, digitlist; base = base)
end

function unfuse_dimensions!(digitlists, digitlist; base::Integer = 2)
    ndim = length(digitlists)
    R = length(digitlist)
    for i = 1:ndim
        for j = 1:R
            digitlists[i][j] =
                _digit_at_index(digitlist[j], ndim - i + 1; numdigits = ndim, base = base)
        end
    end
    return digitlists
end


"""
    interleave_dimensions(digitlists...)

Interleaves the indices of all digitlists into one long digitlist. Use this for
quantics representation of multidimensional objects without fusing indices.
Inverse of [`deinterleave_dimensions`](@ref).
"""
function interleave_dimensions(digitlists...)::Vector{Int}
    results = Vector{Int}(undef, length(digitlists[1]) * length(digitlists))
    return interleave_dimensions!(results, digitlists...)
end

function interleave_dimensions!(
    interleaved_digitlist::AbstractArray{<:Integer},
    digitlists...,
)
    idx = 1
    for i in eachindex(digitlists[1])
        for d in eachindex(digitlists)
            interleaved_digitlist[idx] = digitlists[d][i]
            idx += 1
        end
    end
    return interleaved_digitlist
end

"""
    deinterleave_dimensions(digitlist, d)

Reverses the interleaving of bits, i.e. yields digitlists for each dimension from
a long interleaved digitlist. Inverse of [`interleave_dimensions`](@ref).
"""
function deinterleave_dimensions(digitlist, d)
    return [digitlist[i:d:end] for i = 1:d]
end


function deinterleave_dimensions!(
    deinterleaved_digitlists::AbstractArray{<:AbstractArray{I}},
    digitlist,
) where {I<:Integer}
    d = length(deinterleaved_digitlists)
    for i = 1:d
        @. deinterleaved_digitlists[i] = digitlist[i:d:end]
    end
    return deinterleaved_digitlists
end

"""
    serialize_dimensions(digitlists...)

Sequences the indices of all digitlists into one long digitlist. Use this for
quantics representation of multidimensional objects without fusing indices.
Inverse of [`deserialize_dimensions`](@ref).
"""
function serialize_dimensions(digitlists...)::Vector{Int}
    results = Vector{Int}(undef, length(digitlists[1]) * length(digitlists))
    return serialize_dimensions!(results, digitlists...)
end


function serialize_dimensions!(
    serial_digitlist::AbstractArray{<:Integer},
    digitlists...,
)
    l = length(digitlists[1])
    @inbounds for i in eachindex(digitlists)
        copyto!(@view(serial_digitlist[(i-1)*l+1:i*l]), digitlists[i])
    end
    return serial_digitlist
end

"""
    deserialize_dimensions(digitlist, d)

Reverses the sequencing of bits, i.e. yields digitlists for each dimension from
a long serialized digitlist. Inverse of [`serialize_dimensions`](@ref).
"""
function deserialize_dimensions(digitlist, d)
    l = length(digitlist) ÷ d
    return [digitlist[(i-1)*l+1:i*l] for i = 1:d]
end

function deserialize_dimensions!(
    deserialized_digitlists::AbstractArray{<:AbstractArray{I}},
    digitlist,
) where {I<:Integer}
    d, l = length(deserialized_digitlists), length(deserialized_digitlists[1])
    @inbounds for i in Base.OneTo(d)
        @. deserialized_digitlists[i] = digitlist[(i-1)*l+1:i*l]
    end
    return deserialized_digitlists
end

"""
    meander_dimensions(digitlists...)

Sequences the indices of all digitlists into one long digitlist. Use this for
quantics representation of multidimensional objects without fusing indices.
Inverse of [`demeander_dimensions`](@ref).
"""
function meander_dimensions(digitlists...)::Vector{Int}
    results = Vector{Int}(undef, length(digitlists[1]) * length(digitlists))
    meander_dimensions!(results, digitlists...)
    return results
end


function meander_dimensions!(
    meander_digitlist::AbstractArray{<:Integer},
    digitlists...,
)
    l = length(digitlists[1])
    for (i, digits) in enumerate(digitlists)
        if isodd(i)
            # odd dimension index (1-based): reverse
            @inbounds for j = 1:l
                meander_digitlist[(i-1)*l + j] = digits[l - j + 1]
            end
        else
            # even dimension index: normal order
            copyto!(@view(meander_digitlist[(i-1)*l+1:i*l]), digits)
        end
    end
    return meander_digitlist
end

"""
    demeander_dimensions(digitlist, d)

Reverses the sequencing of bits, i.e. yields digitlists for each dimension from
a long meandered digitlist. Inverse of [`meander_dimensions`](@ref).
"""
function demeander_dimensions(digitlist, d)
    result = Vector{Vector{Int}}(undef, d)
    l = length(digitlist) ÷ d
    @inbounds for i in Base.OneTo(d)
        segment = digitlist[(i-1)*l+1:i*l]
        result[i] = isodd(i) ? reverse(segment) : segment
    end
    return result
end

function demeander_dimensions!(
    demeandered_digitlists::AbstractArray{<:AbstractArray{I}},
    digitlist,
) where {I<:Integer}
    d, l = length(demeandered_digitlists), length(demeandered_digitlists[1])
    @inbounds for i in Base.OneTo(d)
        segment = @view(digitlist[(i-1)*l+1:i*l])
        if isodd(i)
            demeandered_digitlists[i] = reverse(segment)
        else
            demeandered_digitlists[i] = Vector(segment)
        end
    end
    return demeandered_digitlists
end

"""
Convert a fused quantics representation to an unfused quantics interleaved representation
"""
function fused_to_interleaved(
    digitlist::AbstractVector{T},
    d;
    base = 2,
)::Vector{T} where {T<:Integer}
    return interleave_dimensions(unfuse_dimensions(digitlist, d; base = base)...)
end
"""
Convert a fused quantics representation to an unfused quantics serial representation
"""
function fused_to_serial(
    digitlist::AbstractVector{T},
    d;
    base = 2,
)::Vector{T} where {T<:Integer}
    return serialize_dimensions(unfuse_dimensions(digitlist, d; base = base)...)
end
"""
Convert a fused quantics representation to an unfused quantics meander representation
"""
function fused_to_meander(
    digitlist::AbstractVector{T},
    d;
    base = 2,
)::Vector{T} where {T<:Integer}
    return meander_dimensions(unfuse_dimensions(digitlist, d; base = base)...)
end

"""
Convert an unfused quantics interleaved representation to an fused quantics representation
"""
function interleaved_to_fused(
    digitlist::AbstractVector{T};
    base::Integer = 2,
    d::Int = 1,
)::Vector{T} where {T<:Integer}
    fused_digitlist = Vector{T}(undef, length(digitlist) ÷ d)
    fuse_dimensions!(fused_digitlist, deinterleave_dimensions(digitlist, d)...)
    return fused_digitlist
end
"""
Convert an unfused quantics serial representation to an fused quantics representation
"""
function serial_to_fused(
    digitlist::AbstractVector{T};
    base::Integer = 2,
    d::Int = 1,
)::Vector{T} where {T<:Integer}
    fused_digitlist = Vector{T}(undef, length(digitlist) ÷ d)
    fuse_dimensions!(fused_digitlist, deserialize_dimensions(digitlist, d)...)
    return fused_digitlist
end
"""
Convert an unfused quantics meander representation to an fused quantics representation
"""
function meander_to_fused(
    digitlist::AbstractVector{T};
    base::Integer = 2,
    d::Int = 1,
)::Vector{T} where {T<:Integer}
    fused_digitlist = Vector{T}(undef, length(digitlist) ÷ d)
    fuse_dimensions!(fused_digitlist, demeander_dimensions(digitlist, d)...)
    return fused_digitlist
end

"""
    function quantics_to_index_fused(
        digitlist::AbstractVector{<:Integer};
        base::Val{B}=Val(2), dims::Val{d}=Val(1)
    )::NTuple{d,Int} where {B, d}

Convert a d-dimensional index from fused quantics representation to d Integers.

* `base`           base for quantics (default: 2)
* `d`           number of dimensions
* `digitlist`     base-b representation

See also [`quantics_to_index`](https://tensors4fields.gitlab.io/quanticstci.jl/dev/#QuanticsTCI.quantics_to_index).
"""
function quantics_to_index_fused(
    digitlist::AbstractVector{<:Integer};
    base::Integer = 2,
    dims::Val{d} = Val(1),
)::NTuple{d,Int} where {d}
    R = length(digitlist)
    result = ones(MVector{d,Int})

    maximum(digitlist) <= base^d || error("maximum(digitlist) <= base^d")
    minimum(digitlist) >= 0 || error("minimum(digitlist) >= 0")

    for n = 1:R # from the least to most significant digit
        scale = base^(n - 1) # length scale
        tmp = digitlist[R-n+1] - 1
        for i = 1:d # in the order of 1st dim, 2nd dim, ...
            div_, rem_ = divrem(tmp, base)
            result[i] += rem_ * scale
            tmp = div_
        end
    end

    return tuple(result...)
end

"""
    function quantics_to_index(
        digitlist::AbstractVector{<:Integer};
        base::Val{B}=Val(2), dims::Val{d}=Val(1),
        unfoldingscheme::Symbol=:fused
    )::NTuple{d,Int} where {B, d}

Convert a d-dimensional index from fused quantics representation to d Integers.

* `base`           base for quantics (default: 2)
* `d`           number of dimensions
* `digitlist`     base-b representation
* `unfoldingscheme`    Unfolding scheme (fused or interleaved)
"""
function quantics_to_index(
    digitlist::AbstractVector{<:Integer};
    base::Integer = 2,
    dims::Val{d} = Val(1),
    unfoldingscheme::Symbol = :fused,
)::NTuple{d,Int} where {d}
    if unfoldingscheme === :fused
        return quantics_to_index_fused(digitlist, base = base, dims = dims)
    elseif unfoldingscheme === :serial
        return quantics_to_index_fused(
            serial_to_fused(digitlist; base = base, d = d),
            base = base,
            dims = dims,
        )
    elseif unfoldingscheme === :meander
        return quantics_to_index_fused(
            meander_to_fused(digitlist; base = base, d = d),
            base = base,
            dims = dims,
        )
    else
        return quantics_to_index_fused(
            interleaved_to_fused(digitlist; base = base, d = d),
            base = base,
            dims = dims,
        )
    end
end


"""
    function index_to_quantics!(digitlist, index::Integer; base::Integer=2)

* `digitlist`     base-b representation (1d vector)
* `base`           base for quantics (default: 2)
"""
function index_to_quantics!(digitlist, index::Integer; base::Integer = 2)
    numdigits = length(digitlist)
    for i = 1:numdigits
        digitlist[i] = mod(index - 1, base^(numdigits - i + 1)) ÷ base^(numdigits - i) + 1
    end
    return digitlist
end

"""
    index_to_quantics(index::Integer; numdigits=8, base::Integer=2)

Does the same as [`index_to_quantics!`](@ref) but returns a new vector.
"""
function index_to_quantics(index::Integer; numdigits = 8, base::Integer = 2)
    digitlist = Vector{Int}(undef, numdigits)
    return index_to_quantics!(digitlist, index; base = base)
end

"""
    index_to_quantics_fused!(digitlist::AbstractVector{<:Integer}, index::NTuple{D,<:Integer}; base::Integer=2) where {D}

Does the opposite of [`quantics_to_index_fused`](@ref)
"""
function index_to_quantics_fused!(
    digitlist::AbstractVector{<:Integer},
    index::NTuple{D,<:Integer};
    base::Integer = 2,
) where {D}
    R = length(digitlist)
    ndims = length(index)
    digitlist .= 1
    for dim = 1:ndims
        for i = 1:R # from the left to right
            digitlist[i] +=
                (base^(dim - 1)) *
                (_digit_at_index(index[dim], i; numdigits = R, base = base) - 1)
        end
    end
    return digitlist
end

function index_to_quantics_fused(
    index::NTuple{D,<:Integer};
    numdigits::Integer = 8,
    base::Integer = 2,
) where {D}
    digitlist = Vector{Int}(undef, numdigits)
    return index_to_quantics_fused!(digitlist, index; base = base)
end


"""
`index`, `position` and the result are one-based.

`base`: base
`index`: The integer to look at.
`position=1`: Specify the position of the digit to look at. 1 is the most significant (most left) digit.
`numdigits=8`: Specify the number of digits in the number `index`.
"""
function _digit_at_index(index, position; numdigits = 8, base::Integer = 2)
    p_ = numdigits - position + 1
    return mod((index - 1), base^p_) ÷ base^(p_ - 1) + 1
end
