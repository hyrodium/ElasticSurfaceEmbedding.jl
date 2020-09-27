push!(LOAD_PATH, "Modules")
using IntervalSets
using Printf
using BSpline
using ElasticSurfaceEmbedding

@ParametricMapping 𝒑₍₀₎(u)=[cos(u[2])*cosh(u[1]),sin(u[2])*cosh(u[1]),u[1]]

n=9

D=(-π/2..π/2,-π/(4n)..π/(4n))
Settings("Catenoid",up=3,down=-3,right=3,left=-3,mesh=(2n,1),unit=150,slack=true)
InitialConfiguration(D)
NewtonMethodIteration(nip=25)
NewtonMethodIteration(nip=25)
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)
p_Refinement([0,1])
h_Refinement([Float64[],[0.0]])
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)

FinalOutput(unitlength=(50,"mm"))

@ParametricMapping 𝒑₍₀₎(u)=[cos(u[1])*cosh(u[2]),sin(u[1])*cosh(u[2]),u[2]]

id=1
D=(-π..π,(id-1)*π/(2n)..id*π/(2n)) #横方向
Settings("Catenoid-"*(@sprintf "%02d" id),up=1,down=-5,right=3,left=-3,mesh=(4n,1),unit=150,slack=true)
InitialConfiguration(D,n₁=15)
NewtonMethodIteration(nip=25)
NewtonMethodIteration(nip=45)
p_Refinement([0,1])
h_Refinement([Float64[],[0.0]])
NewtonMethodIteration(nip=45)

FinalOutput(unitlength=(50,"mm"))
