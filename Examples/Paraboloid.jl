## Load packages
using IntervalSets
using BasicBSpline
using ElasticSurfaceEmbedding

## Set Parametric mapping
@parametric_mapping 𝒑₍₀₎(u) = [u...,u'*u]
D(i,n) = (-1.0..1.0, (i-1)/n..i/n)
name = "Paraboloid-a"
settings(name,up=2,down=-2,right=2,left=-2,mesh=(20,1),unit=200,slack=true,colorbarsize=0.3)

## Numrcical Computation
i=1
print_strain(D(i,10))
initial_configulation(D(i,10), n₁=19)
initial_configulation(D(i,10))

newton_onestep(fixingmethod=:FixThreePoints)
newton_onestep()
spline_refinement(p₊=[0,1],k₊=[Knots([]),Knots([(i-1/2)/10])])
newton_onestep()
pin_state(tag="$(name)-"*string(i+2))

## Export pinned states
computed_shapes()
print_knots()
export_pinned_states(unitlength=(100,"mm"))
