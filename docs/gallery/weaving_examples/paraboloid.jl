# ---
# title: Paraboloid
# cover: ../assets/paraboloid.jpg
# description: Weaving a simple curved surface with positive curvature.
# ---

# Weaving a simple curved surface with positive curvature.

# ![](../assets/paraboloid.jpg)

# ## Load packages
using IntervalSets
using BasicBSpline
using StaticArrays
using ElasticSurfaceEmbedding

# ## Define the shape of the surface
ElasticSurfaceEmbedding.𝒑₍₀₎(u¹,u²) = SVector(u¹, u², u¹^2+u²^2)
n = 10
D(i,n) = (-1.0..1.0, (i-1)/n..i/n)

# ## Compute the shape of the embeddings
steptree = ElasticSurfaceEmbedding.StepTree()
for i in 1:10
    initial_state!(steptree, D(i,n), n₁=33)
    newton_onestep!(steptree, fixingmethod=:fix3points)
    newton_onestep!(steptree)
    newton_onestep!(steptree)
    newton_onestep!(steptree)
    newton_onestep!(steptree)
    refinement!(steptree, p₊=(0,1), k₊=(EmptyKnotVector(), KnotVector([i/n-1/2n])))
    newton_onestep!(steptree)
    newton_onestep!(steptree)
    pin!(steptree)
end

# ## Export the shape in SVG format
export_pinned_steps("paraboloid", steptree, xlims=(-2,2), ylims=(-2,2), unitlength=(100,"mm"), mesh=(20,1))

# ![](paraboloid/pinned/pinned-9.svg)
# ![](paraboloid/pinned/pinned-18.svg)
# ![](paraboloid/pinned/pinned-27.svg)
# ![](paraboloid/pinned/pinned-36.svg)
# ![](paraboloid/pinned/pinned-45.svg)
# ![](paraboloid/pinned/pinned-54.svg)
# ![](paraboloid/pinned/pinned-63.svg)
# ![](paraboloid/pinned/pinned-72.svg)
# ![](paraboloid/pinned/pinned-81.svg)
# ![](paraboloid/pinned/pinned-90.svg)
