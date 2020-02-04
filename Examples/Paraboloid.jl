using Distributed
addprocs(13);
@everywhere push!(LOAD_PATH, "Modules")
@everywhere push!(LOAD_PATH, ".")
using IntervalSets
using Printf
using Bspline
using I4SM

@DefineShape 𝒑₍₀₎(u)=[u...,u'*(u.*[-1.0,1.0])]
D=((-1.0)..1.0,0.0..0.1)
Settings("Paraboloid-1a",up=3,down=-3,right=3,left=-3,mesh=(20,1),unit=200,slack=true)

InitialConfiguration(D)

fixed(n₁,n₂)=[[1,(n₂+1)÷2,1],[1,(n₂+1)÷2,2],[n₁,(n₂+1)÷2,1],[n₁,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2,1],[(n₁+1)÷2,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2-1,1]]
NewtonMethodIteration(nip=25,fixed=fixed)
NewtonMethodIteration(nip=45)
p_Refinement([0,1])
h_Refinement([Float64[],[0.05]])
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)
