using Distributed
addprocs(1);
@everywhere push!(LOAD_PATH, "Modules")
using IntervalSets
using Printf
using BSpline
using ElasticSurfaceEmbedding

@DefineShape 𝒑₍₀₎(u)=[cos(u[2])*sinh(u[1]),sin(u[2])*sinh(u[1]),u[2]]

D=(-π/2..π/2,-π/(4n)..π/(4n))
Settings("Helicoid",up=3,down=-3,right=3,left=-3,mesh=(2n,1),unit=150,slack=true)
InitialConfiguration(D)
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
Refinement(p₊=[0,1],k₊=[Knots([]),Knots([0.0])])
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
FinalOutput(unitlength=(50,"mm"))


@DefineShape 𝒑₍₀₎(u)=[cos(u[1])*sinh(u[2]),sin(u[1])*sinh(u[2]),u[1]]

n=9
id=7
D=(-π..π,(id-1)*π/(2n)..id*π/(2n)) #横方向
Settings("Helicoid-"*(@sprintf "%02d" id),up=1,down=-5,right=3,left=-3,mesh=(4n,1),unit=150,slack=true)
InitialConfiguration(D,n₁=15)
NewtonMethodIteration()
NewtonMethodIteration()
Refinement(p₊=[0,1],k₊=[Knots([]),Knots([(id-1/2)*π/(2n)])])
NewtonMethodIteration()

FinalOutput(unitlength=(50,"mm"))
