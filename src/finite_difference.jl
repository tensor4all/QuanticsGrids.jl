# Periodic central finite differences on a discretized grid.
#
# Rewritten from the periodic central-difference operators (QTT/MPO form) of
# sakurairihito's QuanticsTCI.jl#9, but expressed on the grid: the spacing is
# taken from `grid_step(g)` instead of being assumed to be `2^-R`, so grids
# whose domain is not the unit interval get the correct coefficient.

"""
    central_difference_periodic(g::DiscretizedGrid{1}; order=1) -> Matrix{Float64}

Return the `n×n` matrix of the periodic central finite difference of the
requested `order` (`1` or `2`) on the one-dimensional grid `g`, where `n` is the
number of grid points of `g`.

Writing `Δ = grid_step(g)` for the grid spacing, the matrix acts on the values
`f_i = f(x_i)` on the grid points `x_i` as

```math
(D_1 f)_i = \\frac{f_{i+1} - f_{i-1}}{2Δ},
\\qquad
(D_2 f)_i = \\frac{f_{i+1} - 2f_i + f_{i-1}}{Δ^2},
```

where the indices are wrapped periodically (`f_0 = f_n`, `f_{n+1} = f_1`).
Because `Δ` is read from the grid, the operator is correct on a grid of any
extent, not only on `[0, 1)`.

The matrix is returned as a dense `Matrix{Float64}`.

# Examples
```julia-repl
julia> g = DiscretizedGrid{1}(2, -2.0, 2.0);   # 4 grid points, spacing Δ = 1

julia> QuanticsGrids.grid_step(g)
1.0

julia> central_difference_periodic(g)
4×4 Matrix{Float64}:
  0.0   0.5   0.0  -0.5
 -0.5   0.0   0.5   0.0
  0.0  -0.5   0.0   0.5
  0.5   0.0  -0.5   0.0

julia> central_difference_periodic(g; order=2)
4×4 Matrix{Float64}:
 -2.0   1.0   0.0   1.0
  1.0  -2.0   1.0   0.0
  0.0   1.0  -2.0   1.0
  1.0   0.0   1.0  -2.0
```

A non-periodic function (for example a linear one) sees a jump in the two rows
that wrap around.

!!! note
    `g` must not include the endpoint (`includeendpoint = false`), since a
    periodic grid identifies its two ends.
"""
function central_difference_periodic(g::DiscretizedGrid{1}; order::Integer=1)
    if !(order in (1, 2))
        throw(ArgumentError(lazy"Got order = $order. Supported orders are 1 and 2."))
    end
    if g.includeendpoint[1]
        throw(ArgumentError(lazy"Got a grid with includeendpoint = true. A periodic grid identifies its two ends, so includeendpoint must be false."))
    end

    n = grid_bases(g)[1]^grid_Rs(g)[1]
    Δ = grid_step(g)
    c = order == 1 ? 1 / (2Δ) : 1 / Δ^2
    D = zeros(Float64, n, n)
    for i in 1:n
        left, right = mod1(i - 1, n), mod1(i + 1, n)
        if order == 1
            D[i, left] -= c
            D[i, right] += c
        else
            D[i, left] += c
            D[i, i] -= 2c
            D[i, right] += c
        end
    end
    return D
end
