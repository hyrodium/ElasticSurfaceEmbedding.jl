using IntervalSets
using BasicBSpline
using StaticArrays
using ElasticSurfaceEmbedding

ElasticSurfaceEmbedding.𝒑₍₀₎(u¹,u²) = SVector(u¹, u², u¹^2+u²^2)
n = 10
D(i,n) = (-1.0..1.0, (i-1)/n..i/n)

steptree = ElasticSurfaceEmbedding.StepTree()
for i in 1:10
    initial_state!(steptree, D(i,n))
    newton_onestep!(steptree, fixingmethod=:fix3points)
    newton_onestep!(steptree)
    newton_onestep!(steptree)
    newton_onestep!(steptree)
    newton_onestep!(steptree)
    refinement!(steptree, p₊=(0,1), k₊=suggest_knotvector(steptree))
    newton_onestep!(steptree)
    newton_onestep!(steptree)
    pin!(steptree)
end

export_pinned_steps("paraboloid", steptree, xlims=(-2,2), ylims=(-2,2), unitlength=(100,"mm"), mesh=(20,1))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
