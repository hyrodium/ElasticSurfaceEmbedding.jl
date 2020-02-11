using Distributed
addprocs(1);
@everywhere push!(LOAD_PATH, "Modules")
using IntervalSets
using Printf
using Bspline
using I4SM

@DefineShape 𝒑₍₀₎(u)=[u...,u'*(u.*[-1.0,1.0])]
n=10
id=1

D=(-1.0..1.0, (i-1)/n..i/n)
Settings("Paraboloid-"*(@sprintf "%02d" id), up=3,down=-3,right=3,left=-3,mesh=(20,1),unit=200,slack=true)
InitialConfiguration(D)
fixed(n₁,n₂)=[[1,(n₂+1)÷2,1],[1,(n₂+1)÷2,2],[n₁,(n₂+1)÷2,1],[n₁,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2,1],[(n₁+1)÷2,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2-1,1]]
NewtonMethodIteration(nip=25,fixed=fixed)
NewtonMethodIteration(nip=45)
p_Refinement([0,1])
h_Refinement([Float64[],[(i-1/2)/n]])
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)

Settings("Paraboloid-10b",up=3,down=-3,right=3,left=-3,mesh=(20,1),unit=200,slack=true)
Restoration()
FinalOutput(unitlength=(50,"mm"),cutout=(0.2,5))
