module QuanticsGrids

export DiscretizedGrid, InherentDiscreteGrid
export quantics_to_grididx, quantics_to_origcoord
export grididx_to_quantics, grididx_to_origcoord
export origcoord_to_quantics, origcoord_to_grididx

abstract type Grid{D} end

include("grid.jl")
include("grid_discretized.jl")

end
