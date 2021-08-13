module ElasticSurfaceEmbedding

using LinearAlgebra
using Printf
using IntervalSets
using JSON

using BasicBSpline
using ExportNURBS

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
