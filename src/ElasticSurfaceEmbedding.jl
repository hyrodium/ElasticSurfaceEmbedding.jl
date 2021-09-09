module ElasticSurfaceEmbedding

using LinearAlgebra
using Printf
using IntervalSets
using JSON

using BasicBSpline
using ExportNURBS

export @parametric_mapping
export settings
export initial_configulation
export newton_onestep
export spline_refinement
export add_pin
export remove_pin
export export_all_pinned_states
export print_strain
export print_knots
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
