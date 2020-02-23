using Distributed
addprocs(1);
@everywhere push!(LOAD_PATH, "Modules")
using IntervalSets
using Printf
using BSpline
using ElasticSurfaceEmbedding

n=30
id=19
@DefineShape 𝒑₍₀₎(u)=5.0*[((1+sqrt(2))*sinpi((u[1]+u[2])))/(-2-sqrt(2)+(1+sqrt(2))*sinpi((u[1]-u[2]))),-(((1+sqrt(2))*cospi((u[1]+u[2])))/(-2-sqrt(2)+(1+sqrt(2))*sinpi((u[1]-u[2])))),(2*(1+sqrt(2))*cot(1/4*π*(3-2u[1]+2u[2])))/(1+(3+2*sqrt(2))*cot(1/4*π*(3-2u[1]+2u[2]))^2)]
D=(0.0..1.0,(id-1)/nf..id/nf)
Settings("Torus"*(@sprintf "%02d" id),up=25,down=-15,right=25,left=-15,mesh=(nf,1),unit=20,slack=true)
InitialConfiguration(D,n₁=33)

fixed(n₁,n₂)=[[1,(n₂+1)÷2,1],[1,(n₂+1)÷2,2],[n₁,(n₂+1)÷2,1],[n₁,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2,1],[(n₁+1)÷2,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2-1,1]]
NewtonMethodIteration(fixed=fixed)
NewtonMethodIteration(fixed=fixed)
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
Refinement(p₊=[0,1],k₊=[Knots([]),Knots([(id-1/2)/n])])
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
