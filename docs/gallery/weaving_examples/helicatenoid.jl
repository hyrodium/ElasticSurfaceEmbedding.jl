# ---
# title: Helicatenoid
# cover: ../assets/helicatenoid.jpg
# description: Weaving a transformable curved surface from catenoid to helicoid.
# ---

# Weaving a transformable curved surface from catenoid to helicoid.

# ```@raw html
# <div class="videoWrapper">
#   <!-- Copy & Pasted from YouTube -->
#   <iframe width="560" height="315" src="https://www.youtube.com/embed/Gp6XkPLCw7s" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>
# </div>
# ```

# ## Load packages
using IntervalSets
using BasicBSpline
using StaticArrays
using ElasticSurfaceEmbedding

# ## Define the shape of the surface (non-periodic direction)
ElasticSurfaceEmbedding.𝒑₍₀₎(u¹,u²) = SVector(cos(u²)*cosh(u¹),sin(u²)*cosh(u¹),u¹)
n = 9
Da(n) = (-π/2..π/2,-π/(4n)..π/(4n))

# ## Compute the shape of the embeddings
show_strain(Da(n))

steptree = initial_state(Da(n), n₁=33)
newton_onestep!(steptree, fixingmethod=:fix3points)
newton_onestep!(steptree)
newton_onestep!(steptree)
newton_onestep!(steptree)
newton_onestep!(steptree)
refinement!(steptree, p₊=(0,1), k₊=(EmptyKnotVector(), KnotVector([0])))
newton_onestep!(steptree)
newton_onestep!(steptree)
pin!(steptree)

# ## Export the shape in SVG format
export_pinned_steps("helicatenoid-a", steptree, unitlength=(40,"mm"), mesh=(18,1))

# ![](helicatenoid-a/pinned/pinned-9.svg)

# ## Define the shape of the surface (periodic direction)
ElasticSurfaceEmbedding.𝒑₍₀₎(u¹,u²) = SVector(cos(u¹)*cosh(u²),sin(u¹)*cosh(u²),u²)
Db(i,n) = (-π..π,(i-1)*π/(2n)..(i)*π/(2n))

## Check the maximum strain
for i in 1:n
    show_strain(Db(i,n))
end

## Numerical computing
steptree = StepTree()
for i in 1:n
    initial_state!(steptree, Db(i,n), n₁=33)
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

# ## Export the shapes in SVG format
export_pinned_steps("helicatenoid-b", steptree, unitlength=(40,"mm"), mesh=(36,1))

# ![](helicatenoid-b/pinned/pinned-9.svg) ![](helicatenoid-b/pinned/pinned-18.svg) ![](helicatenoid-b/pinned/pinned-27.svg)
# ![](helicatenoid-b/pinned/pinned-36.svg) ![](helicatenoid-b/pinned/pinned-45.svg) ![](helicatenoid-b/pinned/pinned-54.svg)
# ![](helicatenoid-b/pinned/pinned-63.svg) ![](helicatenoid-b/pinned/pinned-72.svg) ![](helicatenoid-b/pinned/pinned-81.svg)
