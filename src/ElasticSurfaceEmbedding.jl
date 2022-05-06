module ElasticSurfaceEmbedding

using LinearAlgebra
using Printf
using Dates
import Statistics.mean

using IntervalSets
using StaticArrays
using JSON
import FileIO.load
import FileIO.save
using OffsetArrays
using ForwardDiff
using FastGaussQuadrature
using Colors
using Luxor

using BasicBSpline
using BasicBSplineExporter

# config
export config_dir, config_slack
export settings
# Numerical computing
export initial_state
export newton_onestep!
export refinement!
# Pin related
export pin, unpin, export_all_pinned_states
# utilities
export show_strain, show_knotvector
export computed_shapes

include("_constants.jl")
include("_integral.jl")
include("_graphics.jl")
include("_slack.jl")
include("_bspline.jl")
include("_geometry.jl")
include("_elasticity.jl")
include("_io.jl")
include("_initialstates.jl")
include("_newton.jl")
include("_pin.jl")

end # module
