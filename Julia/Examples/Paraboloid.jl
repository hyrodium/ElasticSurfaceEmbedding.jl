## if use Distibuted
using Distributed
addprocs(10);
@everywhere push!(LOAD_PATH, "Julia/Modules")

## if do not use Distibuted
push!(LOAD_PATH, "Julia/Modules")

##
using Revise
using IntervalSets
using Printf
using BSpline
using ElasticSurfaceEmbedding

##
@ParametricMapping 𝒑₍₀₎(u)=[u...,u'*u]
D(i,n)=(-1.0..1.0, (i-1)/n..i/n)

##
Settings("Paraboloid_p",up=2,down=-2,right=2,left=-2,mesh=(20,1),unit=200,slack=true,colorbarsize=0.3)

i=1
ShowMaximumStrain(D(i,10))
InitialConfiguration(D(i,10))
NewtonMethodIteration(fixingmethod=:FixThreePoints)
NewtonMethodIteration()
SplineRefinement(p₊=[0,1],k₊=[Knots([]),Knots([(1-1/2)/10])])
NewtonMethodIteration()
PinState(tag="paraboloid-"*string(i))

##
RemovePin(15)
ComputedShapes()
ShowKnots()
ExportPinnedStates(unitlength=(100,"mm"))
