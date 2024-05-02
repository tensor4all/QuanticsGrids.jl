# QuanticsGrids

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tensor4all.github.io/QuanticsGrids.jl/dev)
[![CI](https://github.com/tensor4all/QuanticsGrids.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/tensor4all/QuanticsGrids.jl/actions/workflows/CI.yml)

This module contains utilities for interpolations of functions in the quantics TCI / quantics tensor train (QTT) format.

## Installation

---

Until the module is available via `Pkg`, use the following instructions.

1. Clone the repository to some convenient path
```sh
$ cd /convenient/path
$ git clone git@gitlab.com:QuanticsGrids/QuanticsGrids.jl.git
```
2. In a julia REPL, tell julia where you put the downloaded repository.
```julia
julia> using Pkg; Pkg.dev("convenient/path/QuanticsGrids.jl")
```
3. You should now be able to import the module.
```julia
julia> using QuanticsGrids
```
---

Once the module has been published, the following will install QuanticsGrids.jl:

```julia
julia> using Pkg; Pkg.add("QuanticsGrids.jl")
```

## Definition
We first introduce a $B$-base presentation ($B=2, 3, 4, \cdots$).
To avoid confusing, we will use the 1-based indexing of Julia below.
We represent a positive integer $X~(\ge 1)$ as

$$
X= \sum_{i=1}^R (x_i-1) \times B^{R-i+1} + 1,
$$

where $x_i$ is either 1 or 2 and $R$ is the number of digits.
In this library, the $B$-base representation of $X$ is represented by the vector

$$[x_1, \cdots, x_R].$$

This library supports two unfolding schemes (interleaved and fused representations) for handling multiple variables.
As an example, we consider three variables $X$, $Y$ and $Z$.
Their $B$-base representations are given by

$$[x_1, \cdots, x_R], [y_1, \cdots, y_R], [z_1, \cdots, z_R],$$

respectively.
The interleaved representation of these variables reads

$$[x_1, y_1, z_1, x_2, y_2, z_2, \cdots, x_R, y_R, z_R].$$

The fused representation is given by

$$[\alpha_1, \alpha_2, \cdots, \alpha_R],$$

where

$$
\alpha_i = (x_i-1) + B(y_i-1) + B^2 (z_i-1) + 1
$$

with

$$
1 \le \alpha_i \le B^3.
$$

This convention is consistent with the column major used in Julia: At each digit level $i$, the bit for $x$ runs fastest.
The fused representaion generalizes to any number of variables.


## Usage
This package contains two main functionalities:

1. Low-level functions for converting betwen linear and quantics representations
2. High-level interface for creating a grid

The normal users will use the second high-level interface.
We will describe its usage.

## Discretized grid
`DiscretizedGrid` can be used to discretize a $d$-dimensional space.

## Creating a one-dimensional grid

We can create a one-dimensional grid by discretizing $x$ axis on $[0, 1)$ with $R$ bits as

```julia
import QuanticsGrids as QD
xmin = 0.0
xmax = 1.0
R = 4
grid = QD.DiscretizedGrid{1}(R, xmin, xmax)
````

Here, `DiscretizedGrid` takes one parameter `1`, which denotes the dimension of the grid.
There are six functions for translating between different reprenstations:
`grididx` (1-based linear index), `quantics` and `origcoord` (original coordiate, i.e., $x$).
In `origcoord_to_quantics` and `origcoord_to_grididx`, if `origcoord` is out of the grid, the function returns the closest point in the grid.

Example:
```julia
quantics = fill(1, R)
origcoord = 0.0
grididx = 1
@assert QD.quantics_to_grididx(grid, quantics) == grididx
@assert QD.quantics_to_origcoord(grid, quantics) == origcoord
@assert QD.grididx_to_quantics(grid, grididx) == quantics
@assert QD.grididx_to_origcoord(grid, grididx) == origcoord
@assert QD.origcoord_to_quantics(grid, origcoord) == quantics
@assert QD.origcoord_to_grididx(grid, origcoord) == grididx

quantics = fill(2, R)
origcoord = 1-1/2^R
grididx = 2^R
@assert QD.quantics_to_grididx(grid, quantics) == grididx
@assert QD.quantics_to_origcoord(grid, quantics) == origcoord
@assert QD.grididx_to_quantics(grid, grididx) == quantics
@assert QD.grididx_to_origcoord(grid, grididx) == origcoord
@assert QD.origcoord_to_quantics(grid, origcoord) == quantics
@assert QD.origcoord_to_grididx(grid, origcoord) == grididx
```

Optionally, one can include the end point `grid_max` in a grid as

```julia
import QuanticsGrids as QD
xmin = 0.0
xmax = 1.0
R = 4
grid = QD.DiscretizedGrid{1}(R, xmin, xmax; includeendpoint=true)

@assert QD.grididx_to_origcoord(grid, 1) == xmin
@assert QD.grididx_to_origcoord(grid, 2^R) == xmax
````

## Creating a $d$-dimensional grid
A $d$-dimensional grid, where each axis is discretized with $R$ bits, can be generated in a similar way as follows.
As an option, you can choose the fused representation (`:fused`) or the interleaved representation (`:interleaved`).

### fused representation
```julia
import QuanticsGrids as QD
xmin, xmax = 0.0, 1.0
ymin, ymax = 0.0, 1.0
zmin, zmax = 0.0, 1.0
R = 4
grid = QD.DiscretizedGrid{3}(R, (xmin,ymin,zmin), (xmax,ymax,zmax); unfoldingscheme=:fused)

quantics = fill(1, R)
origcoord = (0.0, 0.0, 0.0)
grididx = (1, 1, 1)
@assert QD.quantics_to_grididx(grid, quantics) == grididx
@assert QD.quantics_to_origcoord(grid, quantics) == origcoord
@assert QD.grididx_to_quantics(grid, grididx) == quantics
@assert QD.grididx_to_origcoord(grid, grididx) == origcoord
@assert QD.origcoord_to_quantics(grid, origcoord) == quantics
@assert QD.origcoord_to_grididx(grid, origcoord) == grididx

# Incrementing the least significant fused bits increments the $x$ index.
quantics = vcat(fill(1, R-1), 2) # [1, 1, ..., 1, 2]
origcoord = (1/2^R, 0.0, 0.0)
grididx = (2, 1, 1)
@assert QD.quantics_to_grididx(grid, quantics) == grididx
@assert QD.quantics_to_origcoord(grid, quantics) == origcoord
@assert QD.grididx_to_quantics(grid, grididx) == quantics
@assert QD.grididx_to_origcoord(grid, grididx) == origcoord
@assert QD.origcoord_to_quantics(grid, origcoord) == quantics
@assert QD.origcoord_to_grididx(grid, origcoord) == grididx
```

### Interleaved representation
```julia
import QuanticsGrids as QD
xmin, xmax = 0.0, 1.0
ymin, ymax = 0.0, 1.0
zmin, zmax = 0.0, 1.0
R = 4
grid = QD.DiscretizedGrid{3}(R, (xmin,ymin,zmin), (xmax,ymax,zmax); unfoldingscheme=QD.UnfoldingSchemes.interleaved)

# (x1, y1, z1, ...., xR, yR, zR)
quantics = fill(1, 3R) # length is 3R
origcoord = (0.0, 0.0, 0.0)
grididx = (1, 1, 1)
@assert QD.quantics_to_grididx(grid, quantics) == grididx
@assert QD.quantics_to_origcoord(grid, quantics) == origcoord
@assert QD.grididx_to_quantics(grid, grididx) == quantics
@assert QD.grididx_to_origcoord(grid, grididx) == origcoord
@assert QD.origcoord_to_quantics(grid, origcoord) == quantics
@assert QD.origcoord_to_grididx(grid, origcoord) == grididx

quantics = vcat(fill(1, 3R-3), [2, 1, 1]) # [1, 1, 1, ..., 2, 1, 1]
origcoord = (1/2^R, 0.0, 0.0)
grididx = (2, 1, 1)
@assert QD.quantics_to_grididx(grid, quantics) == grididx
@assert QD.quantics_to_origcoord(grid, quantics) == origcoord
@assert QD.grididx_to_quantics(grid, grididx) == quantics
@assert QD.grididx_to_origcoord(grid, grididx) == origcoord
@assert QD.origcoord_to_quantics(grid, origcoord) == quantics
@assert QD.origcoord_to_grididx(grid, origcoord) == grididx
```

## Inherent discrete grid
`InherentDiscreteGrid` can be used if the target space is discrete.
`InherentDiscreteGrid`  has a very similar interface to `DiscretizedGrid`.
We provide one example.

```julia
import QuanticsGrids as QD
R = 4
# Grid: [0, 1, ..., 2^R-1]. The second argument (0,) specifies the origin.
grid = QD.InherentDiscreteGrid{1}(R, 0; step=1)

quantics = fill(1, R)
origcoord = 0
grididx = 1
@assert QD.quantics_to_grididx(grid, quantics) == grididx
@assert QD.quantics_to_origcoord(grid, quantics) == origcoord
@assert QD.grididx_to_quantics(grid, grididx) == quantics
@assert QD.grididx_to_origcoord(grid, grididx) == origcoord
@assert QD.origcoord_to_quantics(grid, origcoord) == quantics
@assert QD.origcoord_to_grididx(grid, origcoord) == grididx


quantics = fill(2, R)
origcoord = 2^R-1
grididx = 2^R
@assert QD.quantics_to_grididx(grid, quantics) == grididx
@assert QD.quantics_to_origcoord(grid, quantics) == origcoord
@assert QD.grididx_to_quantics(grid, grididx) == quantics
@assert QD.grididx_to_origcoord(grid, grididx) == origcoord
@assert QD.origcoord_to_quantics(grid, origcoord) == quantics
@assert QD.origcoord_to_grididx(grid, origcoord) == grididx
```

## Use a base other than `2`
When creating a grid, we may want to choose a different base other than 2.

```julia
import QuanticsGrids as QD
R = 4
base = 10
# Grid: [0, 1, ..., 10^R-1]
grid = QD.InherentDiscreteGrid{1}(R, (0,); base=10)

quantics = fill(base, R)
origcoord = base^R-1
grididx = base^R
@assert QD.quantics_to_grididx(grid, quantics) == grididx
@assert QD.quantics_to_origcoord(grid, quantics) == origcoord
@assert QD.grididx_to_quantics(grid, grididx) == quantics
@assert QD.grididx_to_origcoord(grid, grididx) == origcoord
@assert QD.origcoord_to_quantics(grid, origcoord) == quantics
@assert QD.origcoord_to_grididx(grid, origcoord) == grididx
```

## Create a function that takes a quantics index as its input
When using `QuanticsGrids.jl` in combination with `TensorCrossInterpolation.jl`,
one can wrap a function to be interpolated to make a fuction that takes a quantics index:

```julia
import QuanticsGrids as QD
import TensorCrossInterpolation as TCI

R = 4
grid = QD.DiscretizedGrid{2}(R, (0.0, 0.0), (1.0, 1.0))

f(x, y) = sin(x + y) # Function to be interpolated

initialpivots = [fill(1, R)]
localdims = fill(2^2, R)
fq = QD.quanticsfunction(Float64, grid, f) # fq takes an quantics index as an input

tci, ranks, errors = TCI.crossinterpolate2(Float64, fq, localdims, initialpivots; tolerance=1e-8)
````

## References
- M. K. Ritter, Y. N. Fernández, M. Wallerberger, J. von Delft, H. Shinaoka, and X. Waintal, *Quantics Tensor Cross Interpolation for High-Resolution, Parsimonious Representations of Multivariate Functions in Physics and Beyond*, [Phys. Rev. Lett. <b>132</b>, 056501 (2024)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.132.056501)/[arXiv:2303.11819](http://arxiv.org/abs/2303.11819).
