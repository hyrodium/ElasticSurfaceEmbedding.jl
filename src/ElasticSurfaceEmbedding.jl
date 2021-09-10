module ElasticSurfaceEmbedding

using LinearAlgebra
using Printf
using Dates
import Statistics.mean

using IntervalSets
using JSON
import FileIO.load
import FileIO.save
using OffsetArrays
using ForwardDiff
using FastGaussQuadrature
using Colors
using Luxor

using BasicBSpline
using ExportNURBS

export @parametric_mapping
export config_dir, config_slack
export settings
export initial_state
export newton_onestep
export spline_refinement
export add_pin
export remove_pin
export export_all_pinned_states
export show_strain
export show_knots
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
