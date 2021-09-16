# [Numerical computation](@id numerical_computation)

## Installation
On Julia's package mode, run the following commands.
```julia-repl
pkg> add IntervalSets
pkg> add BasicBSpline
pkg> add StaticArrays
pkg> add https://github.com/hyrodium/ElasticSurfaceEmbedding.jl
```

## Overview of our method
The computation process proceeds as shown in the following flowchart (from our paper):

![](img/flowchart.png)

## Example: Paraboloid
Through this section, we treat a paraboloid ``z=x^2+y^2`` as an example.

![](img/Paraboloid1.png)

### Load packages, and optional configuration
Load packages with the following script.
```julia
using IntervalSets
using BasicBSpline
using StaticArrays
using ElasticSurfaceEmbedding
```

Set the output directory. *(optional)*
```julia
config_dir("~/ElasticSurfaceEmbedding-Result")
```

```@docs
config_dir
```

Configure the slack bot. *(optional)*
```julia
config_slack(channel="XXXXXXXXXXX",token="xoxb-0000000000000-0000000000000-XXXXXXXXXXXXXXXXXXXXXXXX")
```

!!! info "Slack bot"
    If you would like to use this feature, the bot in the channel must have the following permissions.
    * [`chat.post`](https://api.slack.com/methods/chat.postMessage)
    * [`files.upload`](https://api.slack.com/methods/files.upload)

```@docs
config_slack
```

### Define the shape of surface
```julia
@parametric_mapping 𝒑₍₀₎(u¹,u²) = SVector(u¹, u², u¹^2+u²^2)
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

!!! info "Limitations of a definition of a surface"
    The definition of the parameterized surface must be self-contained.
    This is because, the package will save the definition as string.
    See `"expr"` key in `~/ElasticSurfaceEmbedding-Result/<name>/<name>.json`.

    For example:
    ```julia
    # This is ok
    @parametric_mapping 𝒑₍₀₎(u¹,u²) = SVector(sin(u¹), sin(u²), u¹*u²)

    # This is ng
    f(x,y) = x*y
    @parametric_mapping 𝒑₍₀₎(u¹,u²) = SVector(sin(u¹), sin(u²), f(u¹,u²))
    ```


!!! info "Direction of the surface"
    In the next step, we'll split the surface into elongated strips.
    The domain of each strip should be a rectangular, and the longer direction is `u¹`, and the shorter direction is `u²`.
    The paraboloid has four‐fold symmetry, so we don't have to take care of it.

```@docs
@parametric_mapping
```

### Settings
Before the computation, we need to set the name of the surface, and output graphics region.
```julia
name = "Paraboloid"
settings(name,canvas=(4,4),mesh=(20,1),unit=200,colorbarsize=0.3)
```

```@docs
settings
```

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
    show_strain(D(i,n))
end
```

The output information will be like this:

![](img/example_strain-prediction.png)

!!! tip "Allowable strain"
    Positive number means tension, and negative number means compression.
    Empirically, it is better if the absolute value of strain is smaller than ``0.01``.

```@docs
show_strain
```

### Initial state
If you finished checking the strain prediction, the next step is determine the initial state.

From this section, the computing is done for each piece of the surface.
First, let's calculate for ``i=1``.
```julia
i=1
```

As a first step, let's compute the initial state.
```julia
initial_state(D(i,n), n₁=19)
```

If you've configured a slack bot, you'll get a message like this:

![](img/initial_state.png)

```@docs
initial_state
```

### Newton-Raphson method iteration

```julia
newton_onestep(fixingmethod=:fix3points)
newton_onestep()
```

You can choose the fixing method from below:
* `:default` (default)
* `:fix3points`

```@docs
newton_onestep
```

### Refinement of B-spline manifold

```julia
spline_refinement(p₊=[0,1],k₊=[Knots(),Knots([(i-1/2)/10])])
```

```@docs
spline_refinement
```

```@docs
show_knots
```

### Pin the state
If you finished computing for the strip, it's time to *pin* the state.
This pin📌 will be used for the next final step.

```julia
add_pin(tag="$(name)-$(i)")
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
![](img/catenoid.png)

```julia
## Load packages
using IntervalSets
using BasicBSpline
using StaticArrays
using ElasticSurfaceEmbedding

## Set parametric mapping (x-direction)
@parametric_mapping 𝒑₍₀₎(u¹,u²) = SVector(cos(u²)*cosh(u¹),sin(u²)*cosh(u¹),u¹)
n=9
Dx(n) = (-π/2..π/2,-π/(4n)..π/(4n))

name = "Catenoid-x"
settings(name,canvas=(8,8),mesh=(18,1),unit=200,colorbarsize=0.3)

## Check the maximum strain
show_strain(Dx(n))

## Numerical computing
initial_state(Dx(n), n₁=19)
newton_onestep(fixingmethod=:fix3points)
newton_onestep()
newton_onestep()
spline_refinement(p₊=(0,1),k₊=(Knots(),Knots(0)))
newton_onestep()
newton_onestep()
newton_onestep()
add_pin(tag="$name")
export_all_pinned_states(unitlength=(30,"mm"))


## Set parametric mapping (y-direction)
@parametric_mapping 𝒑₍₀₎(u¹,u²) = SVector(cos(u¹)*cosh(u²),sin(u¹)*cosh(u²),u²)
n=9
Dy(i,n) = (-π..π,(i-1)*π/(2n)..(i)*π/(2n))

name = "Catenoid-y"
settings(name,canvas=(8,8),mesh=(36,1),unit=200,colorbarsize=0.3)

## Check the maximum strain
for i in 1:9
    show_strain(Dy(i,n))
end

## Numerical computing
for i in 1:9
    initial_state(Dy(i,n), n₁=19)
    newton_onestep(fixingmethod=:fix3points)
    newton_onestep()
    newton_onestep()
    spline_refinement(p₊=(0,1),k₊=(Knots(),Knots((i-1/2)*π/(2n))))
    newton_onestep()
    newton_onestep()
    newton_onestep()
    add_pin(tag="$name-$i")
end
export_all_pinned_states(unitlength=(30,"mm"))
```

### Helicoid
![](img/helicoid.png)

```julia
## Load packages
using IntervalSets
using BasicBSpline
using StaticArrays
using ElasticSurfaceEmbedding

## Set parametric mapping (x-direction)
@parametric_mapping 𝒑₍₀₎(u¹,u²) = SVector(cos(u²)*sinh(u¹),sin(u²)*sinh(u¹),u²)
n=9
Dx(n) = (-π/2..π/2,-π/(4n)..π/(4n))

name = "Helicoid-x"
settings(name,canvas=(8,8),mesh=(18,1),unit=200,colorbarsize=0.3)

## Check the maximum strain
show_strain(Dx(n))

## Numerical computing
initial_state(Dx(n), n₁=19)
newton_onestep(fixingmethod=:fix3points)
newton_onestep()
newton_onestep()
spline_refinement(p₊=(0,1),k₊=(Knots(),Knots(0)))
newton_onestep()
newton_onestep()
newton_onestep()
add_pin(tag="$name")
export_all_pinned_states(unitlength=(30,"mm"))


## Set parametric mapping (y-direction)
@parametric_mapping 𝒑₍₀₎(u¹,u²) = SVector(cos(u¹)*sinh(u²),sin(u¹)*sinh(u²),u¹)
n=9
Dy(i,n) = (-π..π,(i-1)*π/(2n)..(i)*π/(2n))

name = "Helicoid-y"
settings(name,canvas=(8,8),mesh=(36,1),unit=200,colorbarsize=0.3)

## Check the maximum strain
for i in 1:9
    show_strain(Dy(i,n))
end

## Numerical computing
for i in 1:9
    initial_state(Dy(i,n), n₁=19)
    newton_onestep(fixingmethod=:fix3points)
    newton_onestep()
    newton_onestep()
    spline_refinement(p₊=(0,1),k₊=(Knots(),Knots((i-1/2)*π/(2n))))
    newton_onestep()
    newton_onestep()
    newton_onestep()
    add_pin(tag="$name-$i")
end
export_all_pinned_states(unitlength=(30,"mm"))
```
