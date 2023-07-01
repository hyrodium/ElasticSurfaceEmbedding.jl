using IntervalSets
using BasicBSpline
using StaticArrays
using ElasticSurfaceEmbedding

ElasticSurfaceEmbedding.𝒑₍₀₎(u¹,u²) = SVector(cos(u²)*cosh(u¹),sin(u²)*cosh(u¹),u¹)
n=9
Dx(n) = (-π/2..π/2,-π/(4n)..π/(4n))

show_strain(Dx(n))

steptree = initial_state(Dx(n), n₁=33)
newton_onestep!(steptree, fixingmethod=:fix3points)
newton_onestep!(steptree)
newton_onestep!(steptree)
newton_onestep!(steptree)
newton_onestep!(steptree)
refinement!(steptree, p₊=(0,1), k₊=(EmptyKnotVector(), KnotVector([0])))
newton_onestep!(steptree)
newton_onestep!(steptree)
pin!(steptree)

export_pinned_steps("helicatenoid-a", steptree, unitlength=(40,"mm"), mesh=(18,1))

ElasticSurfaceEmbedding.𝒑₍₀₎(u¹,u²) = SVector(cos(u¹)*cosh(u²),sin(u¹)*cosh(u²),u²)
Dy(i,n) = (-π..π,(i-1)*π/(2n)..(i)*π/(2n))

# Check the maximum strain
for i in 1:9
    show_strain(Dy(i,n))
end

# Numerical computing
steptree = StepTree()
for i in 1:9
    initial_state!(steptree, Dy(i,n), n₁=33)
    newton_onestep!(steptree, fixingmethod=:fix3points)
    newton_onestep!(steptree)
    newton_onestep!(steptree)
    newton_onestep!(steptree)
    newton_onestep!(steptree)
    refinement!(steptree, p₊=(0,1), k₊=(EmptyKnotVector(),KnotVector([(i-1/2)*π/(2n)])))
    newton_onestep!(steptree)
    newton_onestep!(steptree)
    pin!(steptree)
end

export_pinned_steps("helicatenoid-b", steptree, unitlength=(40,"mm"), mesh=(36,1))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

