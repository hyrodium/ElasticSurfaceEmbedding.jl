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

include("Constants.jl")
include("NumericalIntegral.jl")
include("Slack.jl")
include("TreeStructure.jl")
include("GeometryAndElasticity.jl")
include("InputOutput.jl")
include("InitialConfiguration.jl")
include("NewtonRaphsonMethod.jl")

end # module
