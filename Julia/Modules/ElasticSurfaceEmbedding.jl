module ElasticSurfaceEmbedding

using LinearAlgebra
using Printf
using Distributed
using IntervalSets
using ForwardDiff
using Dates
using DifferentialEquations
using JLD
using JSON

using BSpline
using ParametricDraw

include("_Constants.jl")
include("_NumericalIntegral.jl")
include("_Slack.jl")
include("_TreeStructure.jl")
include("_GeometryAndElasticity.jl")
include("_InputOutput.jl")
include("_InitialConfiguration.jl")
include("_NewtonRaphsonMethod.jl")

end # module
