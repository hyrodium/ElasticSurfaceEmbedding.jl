module ElasticSurfaceEmbedding

using LinearAlgebra
using Printf
using IntervalSets
using JSON
using Distributed

using BasicBSpline

include("_Constants.jl")
include("_NumericalIntegral.jl")
include("_Images.jl")
include("_Slack.jl")
include("_BSpline.jl")
include("_GeometryAndElasticity.jl")
include("_InputOutput.jl")
include("_InitialConfiguration.jl")
include("_NewtonRaphsonMethod.jl")
include("_Pin.jl")

end # module
