module ElasticSurfaceEmbedding

using LinearAlgebra
using Printf
using IntervalSets
using JLD
using Distributed

using BSpline

include("_Constants.jl")
include("_NumericalIntegral.jl")
include("_Slack.jl")
include("_TreeStructure.jl")
include("_GeometryAndElasticity.jl")
include("_InputOutput.jl")
include("_InitialConfiguration.jl")
include("_NewtonRaphsonMethod.jl")

end # module
