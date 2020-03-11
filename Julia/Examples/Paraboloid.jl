# %% if use Distibuted
using Distributed
addprocs(13);
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
D(i,n)=(-1.0..1.0, (i-1)/n..i/n)

# %%
Settings("XX011",up=3,down=-3,right=3,left=-3,mesh=(20,1),unit=200,slack=true)
InitialConfiguration(D(1,10))
NewtonMethodIteration(fixingmethod=:FixThreePoints)
# NewtonMethodIteration(parent=14)
Refinement(p₊=[0,1],k₊=[Knots([]),Knots([(id-1/2)/10])],parent=2)
NewtonMethodIteration(parent=1)
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()
NewtonMethodIteration()

ShowMaximumStrain(D(1,10))

# %%
@ParametricMapping 𝒑₍₀₎(u)=[u...,u[1]*u[2]]
D(i,n)=(-1.0..1.0, (i-1)/n..i/n)

# %%
Settings("XXX012",up=3,down=-3,right=3,left=-3,mesh=(20,1),unit=200,slack=true)

𝒑₍₀₎([1,1])
