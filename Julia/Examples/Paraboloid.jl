## Load packages
using IntervalSets
using BasicBSpline
using ElasticSurfaceEmbedding

## Set Parametric mapping
@ParametricMapping 𝒑₍₀₎(u) = [u...,u'*u]
D(i,n) = (-1.0..1.0, (i-1)/n..i/n)
name = "Paraboloid"
Settings(name,up=2,down=-2,right=2,left=-2,mesh=(20,1),unit=200,slack=true,colorbarsize=0.3)

## Numrcical Computation
i=1
ShowMaximumStrain(D(i,10))
InitialConfiguration(D(i,10), n₁=19)
InitialConfiguration(D(i,10))

NewtonMethodIteration(fixingmethod=:FixThreePoints)
NewtonMethodIteration()
SplineRefinement(p₊=[0,1],k₊=[Knots([]),Knots([(i-1/2)/10])])
NewtonMethodIteration()
PinState(tag="$(name)-"*string(i+2))

## Export pinned states
ComputedShapes()
ShowKnots()
ExportPinnedStates(unitlength=(100,"mm"))
