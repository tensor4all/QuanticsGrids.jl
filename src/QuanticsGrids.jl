module QuanticsGrids

export DiscretizedGrid, InherentDiscreteGrid
export quantics_to_grididx, quantics_to_origcoord
export grididx_to_quantics, grididx_to_origcoord
export origcoord_to_quantics, origcoord_to_grididx

abstract type Grid{D} end

include("inherent_discrete_grid.jl")
include("discretized_grid.jl")
include("show_grids.jl")

end
