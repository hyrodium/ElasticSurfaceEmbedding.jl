using Distributed
addprocs(14);
@everywhere push!(LOAD_PATH, "Modules")
using IntervalSets
using Printf
using Bspline
using I4SM

@DefineShape 𝒑₍₀₎(u)=[u...,u'*(u.*[-1.0,1.0])]
n=10
id=3


D=(-1.0..1.0, (id-1)/10..id/10)
Settings("Paraboloid-"*(@sprintf "%02d" id),up=3,down=-3,right=3,left=-3,mesh=(20,1),unit=200,slack=true)
InitialConfiguration(D)
fixed(n₁,n₂)=[[1,(n₂+1)÷2,1],[1,(n₂+1)÷2,2],[n₁,(n₂+1)÷2,1],[n₁,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2,1],[(n₁+1)÷2,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2-1,1]]
NewtonMethodIteration(fixed=fixed)
NewtonMethodIteration()
p_Refinement([0,1])
h_Refinement([Float64[],[0.05]])
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
