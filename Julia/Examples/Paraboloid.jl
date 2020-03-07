# %% if use Distibuted
using Distributed
addprocs(15);
@everywhere push!(LOAD_PATH, "Julia/Modules")

# %% if do not use Distibuted
push!(LOAD_PATH, "Julia/Modules")

# %%
using Revise
using IntervalSets
using Printf
using BSpline
using ElasticSurfaceEmbedding

# %%
@ParametricMapping 𝒑₍₀₎(u)=[u...,u'*u]
n=10
id=1

# %%
D=(-1.0..1.0, (id-1)/n..id/n)
Settings("PParaboloid-"*(@sprintf "%02d" id),up=3,down=-3,right=3,left=-3,mesh=(20,1),unit=200,slack=true,overwrite=true)
InitialConfiguration(D)
fixed(n₁,n₂)=[[1,(n₂+1)÷2,1],[1,(n₂+1)÷2,2],[n₁,(n₂+1)÷2,1],[n₁,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2,1],[(n₁+1)÷2,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2-1,1]]
NewtonMethodIteration(fixed=fixed)
NewtonMethodIteration()
Refinement(p₊=[0,1],k₊=[Knots([]),Knots([(id-1/2)/10])])
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()

ShowMaximumStrain(D)
