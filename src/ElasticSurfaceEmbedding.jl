module ElasticSurfaceEmbedding

using LinearAlgebra
using Printf
using Dates
import Statistics.mean

using IntervalSets
using StaticArrays
import FileIO.load
import FileIO.save
using OffsetArrays
using ForwardDiff
using FastGaussQuadrature
using Colors
using Luxor
using ColorBlendModes

using BasicBSpline
using BasicBSplineExporter

# Numerical computing
export initial_state, initial_state!
export newton_onestep!
export refinement!
# Pin related
export pin, unpin
# Exports
export export_all_steps, export_pinned_steps
# utilities
export show_strain, show_knotvector

include("_constants.jl")
include("_graphics.jl")
include("_bspline.jl")
include("_geometry.jl")
include("_elasticity.jl")
include("_io.jl")
include("_initialstates.jl")
include("_newton.jl")
include("_pin.jl")

end # module
