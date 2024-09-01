# [Numerical computation](@id numerical_computation)

## Installation
On Julia's package mode, run the following commands.
```julia-repl
pkg> add IntervalSets
pkg> add StaticArrays
pkg> add BasicBSpline
pkg> add BasicBSplineExporter.jl
pkg> add ElasticSurfaceEmbedding.jl
```

## Overview of our method
Our theoretical framework is based on:

* Mathematical model: Nonlinear elasticity on Riemannian manifold
* Geometric representation: B-spline manifold
* Numerical analysis: Galerkin method, Newton-Raphson method

The computation process proceeds as shown in the following flowchart (from our paper):

![](img/flowchart.png)

For more information, read [our paper](https://arxiv.org/abs/2211.06372) or contact [me](https://twitter.com/Hyrodium)!

## Example: Paraboloid
Through this section, we treat a paraboloid ``z=x^2+y^2`` as an example.

![](img/Paraboloid1.png)

### Load packages, and optional configuration
Load packages with the following script.
```@example paraboloid
using IntervalSets
using BasicBSpline
using StaticArrays
using ElasticSurfaceEmbedding
```

### Define the shape of surface
```@example paraboloid
ElasticSurfaceEmbedding.𝒑₍₀₎(u¹,u²) = SVector(u¹, u², u¹^2+u²^2)
```

```math
\begin{aligned}
\bm{p}_{[0]}(u^1, u^2)
&= \begin{pmatrix}
u^1 \\
u^2 \\
(u^1)^2 + (u^2)^2
\end{pmatrix} \\
D
&= [-1,1]\times[-1,1]
\end{aligned}
```

!!! info "Direction of the surface"
    In the next step, we will split the surface into elongated strips.
    The domain of each strip should be rectangular, and the longer direction is `u¹`, and the shorter direction is `u²`.
    The paraboloid has four‐fold symmetry, so we don't have to take care of it.

### Split the surface into strips
The domain ``D`` will be split into ``D_i``.

```math
\begin{aligned}
D_i
&= [-1,1]\times\left[\frac{i-1}{10},\frac{i}{10}\right] & (i=1,\dots,10)
\end{aligned}
```

![](img/Paraboloid2.png)

In julia script, just define a domain of the strip with function `D(i,n)`.

```@example paraboloid
n = 10
D(i,n) = (-1.0..1.0, (i-1)/n..i/n)
```

### Check the strain prediction
Before computing the embedding numerically, we can predict the strain with *Strain Approximation Formula*:

```math
\begin{aligned}
E_{11}^{\langle 0\rangle}&\approx\frac{1}{2}K_{[0]}B^2\left(r^2-\frac{1}{3}\right)
\end{aligned}
```

You can check this strain estimation using the [`show_strain`](@ref) function.

```@example paraboloid
for i in 1:n
    show_strain(D(i,n))
end
```

!!! tip "Allowable strain"
    Positive number means tension, and negative number means compression.
    Empirically, it is better if the absolute value of strain is smaller than ``0.01 (=1\%)``.

### Initial state
If you finished checking the strain prediction, the next step is determination of the initial state with [`initial_state`](@ref) (or [`initial_state!`](@ref) from the second time).

From this section, the computing is done for each piece of the surface.
First, let's calculate for ``i=1``.
```@example paraboloid
i = 1
```

As a first step, let's compute the initial state.

```@example paraboloid
steptree = initial_state(D(i,n))
```

### Newton-Raphson method iteration

[`newton_onestep!`](@ref) function calculates one step of Newton-Raphson method iteration.

```@example paraboloid
newton_onestep!(steptree, fixingmethod=:fix3points)
newton_onestep!(steptree)
```

You can choose the fixing method from below:
* `:default` (default)
* `:fix3points`

### Refinement of B-spline manifold

```@example paraboloid
refinement!(steptree, p₊=(0,1), k₊=(EmptyKnotVector(),KnotVector([(i-1/2)/10])))
```

The knotvector to be inserted in [`refinement!`](@ref) can be suggested by [`show_knotvector`](@ref) function.

### Pin the step
If you finished computing for the strip, it's time to *pin* the step.
This [`pin!`](@ref) function will be used for the the final export step.

```@example paraboloid
pin!(steptree)
```

If you add a pin mistakenly, you can remove the pin with [`unpin!`](@ref) function.

```@example paraboloid
unpin!(steptree, 4)
```

### Compute more
```@example paraboloid
newton_onestep!(steptree)
newton_onestep!(steptree)
pin!(steptree)

i = 2
initial_state!(steptree, D(i,n))
newton_onestep!(steptree, fixingmethod=:fix3points)
newton_onestep!(steptree)
refinement!(steptree, p₊=(0,1), k₊=(EmptyKnotVector(),KnotVector([(i-1/2)/10])))
newton_onestep!(steptree)
newton_onestep!(steptree)
pin!(steptree)
```

### Export all pinned shapes
This is the final step of the computational process with [`export_pinned_steps`](@ref).

```@example paraboloid
export_pinned_steps(".", steptree, unitlength=(50, "mm"), mesh=(20,1), xlims=(-2,2), ylims=(-0.3,0.3))
```

This will create SVG files in `./pinned`.

`pinned/pinned-6.svg`
![](pinned/pinned-6.svg)

`pinned/pinned-12.svg`
![](pinned/pinned-12.svg)

The all outputs for `i in 1:10` will be like this:
![](img/Paraboloid3.png)

You can edit these files, and craft them into curved surface shape.

## Other examples
can be found in [gallery](@ref gallery)
