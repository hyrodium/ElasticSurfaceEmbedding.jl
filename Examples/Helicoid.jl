using IntervalSets
using BSpline
using ElasticSurfaceEmbedding

n=9
domain = (-π/2..π/2,-π/(4n)..π/(4n))
@parametric_mapping 𝒑₍₀₎(u)=[cos(u[2])*sinh(u[1]),sin(u[2])*sinh(u[1]),u[2]]
@parametric_mapping 𝒑₍₀₎(u)=[cos(u[1])*sinh(u[2]),sin(u[1])*sinh(u[2]),u[1]]
