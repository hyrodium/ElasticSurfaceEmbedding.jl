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

D(i) = (-1..1, (i-1)/10..i/10)
name = "CutoffHexagon"
Settings(name,up=5,down=-5,right=8,left=-8,mesh=(20,1),unit=100,slack=false,colorbarsize=0.3)

## Numrcical Computation
i = 2
ShowMaximumStrain(D(i))
InitialConfiguration(D(i), n₁=35)

NewtonMethodIteration(fixingmethod=:FixThreePoints)
NewtonMethodIteration()
NewtonMethodIteration()
SplineRefinement(p₊=[0,1],k₊=[Knots(),Knots([(i-1/2)/10])])
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
PinState(tag="$(name)-"*string(i))
