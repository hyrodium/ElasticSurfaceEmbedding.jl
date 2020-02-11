using Distributed
addprocs(1);
@everywhere push!(LOAD_PATH, "Modules")
using IntervalSets
using Printf
using Bspline
using I4SM

@DefineShape 𝒑₍₀₎(u)=5.0*[((1+sqrt(2))*sinpi((u[1]+u[2])))/(-2-sqrt(2)+(1+sqrt(2))*sinpi((u[1]-u[2]))),-(((1+sqrt(2))*cospi((u[1]+u[2])))/(-2-sqrt(2)+(1+sqrt(2))*sinpi((u[1]-u[2])))),(2*(1+sqrt(2))*cot(1/4*π*(3-2u[1]+2u[2])))/(1+(3+2*sqrt(2))*cot(1/4*π*(3-2u[1]+2u[2]))^2)]
n=30
id=1

D=(0.0..1.0,(id-1)/n..id/n)
Settings("Torus-"*(@sprintf "%02d" id),up=25,down=-15,right=25,left=-15,mesh=(n,1),unit=20,slack=true)
InitialConfiguration(D,n₁=33)
fixed(n₁,n₂)=[[1,(n₂+1)÷2,1],[1,(n₂+1)÷2,2],[n₁,(n₂+1)÷2,1],[n₁,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2,1],[(n₁+1)÷2,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2-1,1]]
NewtonMethodIteration(nip=25,fixed=fixed)
NewtonMethodIteration(nip=25,fixed=fixed)
NewtonMethodIteration(nip=25)
NewtonMethodIteration(nip=25)
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)
p_Refinement([0,1])
h_Refinement([Float64[],[(id-1/2)/nf]])
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)
