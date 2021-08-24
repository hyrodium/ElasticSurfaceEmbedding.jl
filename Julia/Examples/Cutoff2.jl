## Load packages
using IntervalSets
using BasicBSpline
using ElasticSurfaceEmbedding

## Set Parametric mapping
@ParametricMapping function 𝒑₍₀₎(u)
    z = 8*(u[1] + im * u[2])^(2/3)
    x, y = z.re, z.im
    s=(sqrt(x^2+y^2)-7.33/2)/(6.0-7.33/2)
    Cf1=ifelse(s>0.0,exp(-sqrt(3/4)/s),0.0)
    Cf2=ifelse(1-s>0.0,exp(-sqrt(3/4)/(1-s)),0.0)
    Cg=Cf1/(Cf1+Cf2)
    return [x,y,1.6*Cg]
end

𝒑₍₀₎([-1,1])

ElasticSurfaceEmbedding.𝒑₍₀₎([1,1])-ElasticSurfaceEmbedding.𝒑₍₀₎([-1,1])

ElasticSurfaceEmbedding.𝒑₍₀₎([1,0])
ElasticSurfaceEmbedding.𝒑₍₀₎([1,-1])-ElasticSurfaceEmbedding.𝒑₍₀₎([1,1])

Dom(i) = (-1..1, (i-1)/10..i/10)

Settings("Cutoff_2-c",up=5,down=-5,right=8,left=-8,mesh=(20,1),unit=100,slack=false,colorbarsize=0.3)

## Numrcical Computation
i = 2
ShowMaximumStrain(Dom(i))
InitialConfiguration(Dom(i), n₁=35)

NewtonMethodIteration(fixingmethod=:FixThreePoints)
NewtonMethodIteration()
NewtonMethodIteration()
SplineRefinement(p₊=[0,1],k₊=[Knots(),Knots([(i-1/2)/10])])
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
PinState(tag="paraboloid-"*string(i))

PinState(tag="paraboloid-"*string(6), parent=49)


ShowKnots(index=94)


k₁₊_tmp = [-0.984375, -0.953125, -0.921875, -0.890625, -0.859375, -0.828125, -0.796875, -0.765625, -0.734375, -0.703125, -0.671875, -0.640625, -0.609375, -0.578125, -0.546875, -0.515625, -0.484375, -0.453125, -0.421875, -0.390625, -0.359375, -0.328125, -0.296875, -0.265625, -0.234375, -0.203125, -0.171875, -0.140625, -0.109375, -0.078125, -0.046875, -0.015625, 0.015625, 0.046875, 0.078125, 0.109375, 0.140625, 0.171875, 0.203125, 0.234375, 0.265625, 0.296875, 0.328125, 0.359375, 0.390625, 0.421875, 0.453125, 0.484375, 0.515625, 0.546875, 0.578125, 0.609375, 0.640625, 0.671875, 0.703125, 0.734375, 0.765625, 0.796875, 0.828125, 0.859375, 0.890625, 0.921875, 0.953125, 0.984375]
k₁₊ = Knots(k₁₊_tmp[(t->0.2<abs(t)<0.7).(k₁₊_tmp)])
SplineRefinement(p₊=[0,0],k₊=[k₁₊,Knots()], parent=94)








RemovePin(70)
ShowKnots()

k₁₊ = Knots([-0.96875, -0.90625, -0.84375, -0.78125, -0.71875, -0.65625, -0.59375, -0.53125, -0.46875, -0.40625, -0.34375, -0.28125, -0.21875, -0.15625, -0.09375, -0.03125, 0.03125, 0.09375, 0.15625, 0.21875, 0.28125, 0.34375, 0.40625, 0.46875, 0.53125, 0.59375, 0.65625, 0.71875, 0.78125, 0.84375, 0.90625, 0.96875])

SplineRefinement(p₊=[0,1],k₊=[k₁₊,Knots()], parent=59)


PinState(tag="paraboloid-"*string(i))


NewtonMethodIteration()

Settings("Cutoff2i",up=5,down=-5,right=8,left=-8,mesh=(16,1),unit=100,slack=false,colorbarsize=0.3)

## Numrcical Computation
Settings("Cutoff2-b_",up=5,down=-5,right=8,left=-8,mesh=(9,1),unit=100,slack=true,colorbarsize=0.3)
i = 1
Dom(i) = (1/10..1, (i-1)/10..i/10)
ShowMaximumStrain(Dom(i))
InitialConfiguration(Dom(i), n₁=55)

NewtonMethodIteration(fixingmethod=:FixThreePoints)
NewtonMethodIteration()
NewtonMethodIteration()
SplineRefinement(p₊=[0,1],k₊=[Knots(),Knots([(i-1/2)/10])])
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
PinState(tag="paraboloid-"*string(i))

ExportPinnedStates(unitlength=(100,"mm"))

