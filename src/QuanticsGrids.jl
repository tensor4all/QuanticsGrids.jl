module QuanticsGrids

using StaticArrays

export DiscretizedGrid, InherentDiscreteGrid
export quantics_to_grididx, quantics_to_origcoord
export grididx_to_quantics, grididx_to_origcoord
export origcoord_to_quantics, origcoord_to_grididx

include("quantics.jl")
include("grid.jl")
include("newgrid_discrete.jl")
include("newgrid.jl")

end
