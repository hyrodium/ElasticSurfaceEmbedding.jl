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

include("_Constants.jl")
include("_NumericalIntegral.jl")
include("_Images.jl")
include("_Slack.jl")
include("_BSpline.jl")
include("_Geometry.jl")
include("_Elasticity.jl")
include("_InputOutput.jl")
include("_InitialConfiguration.jl")
include("_NewtonRaphsonMethod.jl")
include("_Pin.jl")

end # module
