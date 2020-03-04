using Distributed
addprocs(1);
@everywhere push!(LOAD_PATH, "Modules")
using IntervalSets
using Printf
using BSpline
using ElasticSurfaceEmbedding
@everywhere using LinearAlgebra

@ParametricMapping function 𝒑₍₀₎(u)
    s=(norm(u)-7.33/2)/(9.0-7.33/2)
    Cf1=ifelse(s>0.0,exp(-sqrt(3/4)/s),0.0)
    Cf2=ifelse(1-s>0.0,exp(-sqrt(3/4)/(1-s)),0.0)
    Cg=Cf1/(Cf1+Cf2)
    return [u[1],u[2],1.6*Cg]
end
n=12
id=5

D=(-12.0..12.0, (id-1)..id)
Settings("Cutoff-"*(@sprintf "%02d" id), up=5,down=-5,right=15,left=-15,mesh=(24,1),unit=200,slack=true)
InitialConfiguration(D,n₁=27)
fixed(n₁,n₂)=[[1,(n₂+1)÷2,1],[1,(n₂+1)÷2,2],[n₁,(n₂+1)÷2,1],[n₁,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2,1],[(n₁+1)÷2,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2-1,1]]
NewtonMethodIteration(,fixed=fixed)
NewtonMethodIteration()
Refinement(p₊=[0,1],k₊=[Knots([]),Knots([id-/2])])
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()

ShowKnots()

FinalOutput(unitlength=(10,"mm"),cutout=(0.2,5))
