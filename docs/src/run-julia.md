# [Numerical computation](@id numerical_computation)

## Installation
On Julia's package mode,
```julia-repl
pkg> add IntervalSets
pkg> add BasicBSpline
pkg> add https://github.com/hyrodium/ElasticSurfaceEmbedding.jl
```

## Overview of our method
The computation process proceeds as shown in the following flowchart (from our paper):

![](img/flowchart.png)

1. Define the shape of surface
1. Split the surface into strips
1. Check the strain prediction
1. Initial configuration
1. Newton-Raphson method iteration
1. Refinement B-spline manifold
1. Newton-Raphson method iteration, again
1. Finish Computation

## Example: Paraboloid
Through this section, we treat a paraboloid ``z=x^2+y^2`` as an example.

![](img/Paraboloid1.png)

### Load packages
```julia
using IntervalSets
using BasicBSpline
using ElasticSurfaceEmbedding
```

### Define the shape of surface
```julia
@parametric_mapping 𝒑₍₀₎(u) = [u[1],u[2],u[1]^2+u[2]^2]
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

```@docs
@parametric_mapping
```

### Settings
#### Output directory (optional)

#### Name and graphics region
```julia
settings(name,up=2,down=-2,right=2,left=-2,mesh=(20,1),unit=200,slack=true,colorbarsize=0.3)
```

```@docs
settings
```

#### Slack (optional)

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

```julia
n = 10
D(i,n) = (-1.0..1.0, (i-1)/n..i/n)
```

### Check the strain prediction
Before numerical computation, we can predict the strain with *Strain Approximation Formula*:

```math
\begin{aligned}
E_{11}^{\langle 0\rangle}&\approx\frac{1}{2}K_{[0]}B^2\left(r^2-\frac{1}{3}\right)
\end{aligned}
```

```julia
for i in 1:n
    print_strain(D(i,n))
end
```

The output information will be like this:

![](img/example_strain-prediction.png)

Positive number means tension, and negative number means compression.
Empirically, it is better if the absolute value of strain is smaller than ``0.01``.

```@docs
print_strain
```

### Initial configuration
If you finished checking the strain prediction, the next step is determine the initial configuration.

```julia
initial_configulation(D(i,n), n₁=19)
```

```@docs
initial_configulation
```

### Newton-Raphson method iteration

```julia
newton_onestep(fixingmethod=:FixThreePoints)
newton_onestep()
```

You can choose the fixing method from below:
* `:DefaultOrientation` (default)
* `:FixThreePoints`

```@docs
newton_onestep
```

### Refinement B-spline manifold

```julia
spline_refinement(p₊=[0,1],k₊=[Knots(),Knots([(i-1/2)/10])])
```

```@docs
spline_refinement
```

```@docs
print_knots
```

### Pin the state
If you finished computing for the strip, it's time to *pin* the state.
This pin📌 will be used for the next final step.

```julia
add_pin(tag="paraboloid-"*string(i))
```

```@docs
add_pin
```

If you add a pin mistakenly, you can remove the pin with `remove_pin` function.

```julia
remove_pin(tag="paraboloid-"*string(i))
remove_pin(10)
```

```@docs
remove_pin
```

### Export all pinned shapes
This is the final step of the computational process.

```@docs
export_all_pinned_states
```

This will create SVG files in `~/ElasticSurfaceEmbedding-Result/Paraboloid/pinned/`.

You can edit these files, and craft them into curved surface shape.

![](img/Paraboloid3.png)

### Utilities
If you want to resume the computation, you can just call `settings` like this:

```julia
settings("Paraboloid")
```

To get the list of computed shapes to resume, `computed_shapes` function can be used.

```julia
computed_shapes()
```

```@docs
computed_shapes
```

## Other examples
### Catenoid
```julia
using IntervalSets
using BasicBSpline
using ElasticSurfaceEmbedding

n=9
domain = (-π/2..π/2,-π/(4n)..π/(4n))
@parametric_mapping 𝒑₍₀₎(u)=[cos(u[2])*cosh(u[1]),sin(u[2])*cosh(u[1]),u[1]]
@parametric_mapping 𝒑₍₀₎(u)=[cos(u[1])*cosh(u[2]),sin(u[1])*cosh(u[2]),u[2]]
```

### Helicoid

