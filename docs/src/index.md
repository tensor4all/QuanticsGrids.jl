# QuanticsGrids.jl user guide

```@meta
CurrentModule = QuanticsGrids
```

This module allows easy translation of functions to quantics representation. It meshes well with the `TensorCrossInterpolation.jl` module, together with which it provides quantics TCI functionality.

# Quickstart

The easiest way to construct a quantics tensor train is the `quanticscrossinterpolate` function. For example, the function ``f(x, y) = \exp(-x - 2y)`` can be interpolated as follows.

```@example simple
using QuanticsGrids
f(x, y) = (cos(x) - cos(x - 2y)) * abs(x + y)
xvals = range(-6, 6; length=256)
yvals = range(-12, 12; length=256)
qtt, ranks, errors = quanticscrossinterpolate(Float64, f, [xvals, yvals]; tolerance=1e-8)
```

The QTT can then be evaluated a function of indices enumerating the `xvals` and `yvals` arrays.

```@example simple
using Plots
qttvals = qtt.(1:256, collect(1:256)')
contour(xvals, yvals, qttvals, fill=true)
xlabel!("x")
ylabel!("y")
savefig("simple.svg"); nothing # hide
```

![](simple.svg)

# Quantics representation

## One dimension  / variable
Functions that translate between "normal" and quantics representation, i.e. between ``i`` and ``u_k`` in
```math
i = \frac{u_1}{1} + \frac{u_2}{2} + \ldots + \frac{u_n}{2^{n-1}}
```
This is effectively just a representation in ``n`` bits, where each bit is a tensor leg of the tensor to be represented.

```@docs
QuanticsGrids.index_to_quantics
QuanticsGrids.quantics_to_index
```

## Multiple dimensions / variables
There are two ways to represent ``D``-dimensional functions in quantics:
1. **Fused representation**: The ``D`` tensor legs corresponding to the same length scale, with local dimension ``d`` (usually ``d=2``) are merged to one leg with dimension ``d^D``.
2. **Interleaved representation**: Tensor legs are interleaved, such that the ``D`` tensor legs corresponding to each length scale are neighbours.

### Translation

To translate between multiple indices and their quantics representation in *fused* representation, use the following functions:
```@docs
QuanticsGrids.index_to_quantics_fused
QuanticsGrids.quantics_to_index_fused
```


To translate between multiple indices and their quantics representation in *interleaved* representation, use the following functions:
```@docs
QuanticsGrids.index_to_quantics_interleaved
QuanticsGrids.quantics_to_index_interleaved
```

### Fusing legs
To fuse legs by hand:
```@docs
QuanticsGrids.fuse_dimensions
QuanticsGrids.split_dimensions
```

### Interleaving legs
To interleave legs:
```@docs
QuanticsGrids.interleave_dimensions
QuanticsGrids.deinterleave_dimensions
```