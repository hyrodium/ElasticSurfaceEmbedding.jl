using IntervalSets
using BSpline
using ElasticSurfaceEmbedding

n=9
domain = (-π/2..π/2,-π/(4n)..π/(4n))
@ParametricMapping 𝒑₍₀₎(u)=[cos(u[2])*cosh(u[1]),sin(u[2])*cosh(u[1]),u[1]]
@ParametricMapping 𝒑₍₀₎(u)=[cos(u[1])*cosh(u[2]),sin(u[1])*cosh(u[2]),u[2]]
