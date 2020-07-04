## If use Distibuted
using Distributed
addprocs(14);
@everywhere push!(LOAD_PATH, "Julia/Modules")

## If not use Distibuted
push!(LOAD_PATH, "Julia/Modules")

## Load packages
using Revise
using IntervalSets
using Printf
using BasicBSpline
using ElasticSurfaceEmbedding

## Set Parametric mapping
@ParametricMapping 𝒑₍₀₎(u) = [u...,u'*u]
D(i,n) = (-1.0..1.0, (i-1)/n..i/n)
Settings("Paraboloid_B",up=2,down=-2,right=2,left=-2,mesh=(20,1),unit=200,slack=true,colorbarsize=0.3)

## Numrcical Computation

i=1
ShowMaximumStrain(D(i,10))
InitialConfiguration(D(i,10), n₁=19)
InitialConfiguration(D(i,10))

NewtonMethodIteration(fixingmethod=:FixThreePoints)
NewtonMethodIteration()
SplineRefinement(p₊=[0,1],k₊=[Knots([]),Knots([(i-1/2)/10])])
NewtonMethodIteration()
PinState(tag="paraboloid-"*string(i+1))

## Export pinned states
RemovePin(16)
ComputedShapes()
ShowKnots()
ExportPinnedStates(unitlength=(100,"mm"))
# TODO: use unitful.jl ExportPinnedStates(unitlength=(100,"mm"))
