# Elastic Surface Embedding; Weaving Parer Strips

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hyrodium.github.io/ElasticSurfaceEmbedding.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://hyrodium.github.io/ElasticSurfaceEmbedding.jl/dev)
[![Build Status](https://github.com/hyrodium/ElasticSurfaceEmbedding.jl/workflows/CI/badge.svg)](https://github.com/hyrodium/ElasticSurfaceEmbedding.jl/actions)
[![Coverage](https://codecov.io/gh/hyrodium/ElasticSurfaceEmbedding.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/hyrodium/ElasticSurfaceEmbedding.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![arXiv](https://img.shields.io/badge/math.DG-arXiv%3A2211.06372-B31B1B.svg)](https://arxiv.org/abs/2211.06372)

## JuliaCon2023 Talk! 📣🕙
I gave a lightning talk about this repository at JuliaCon2023!

[![JuliaCon2023 Talk](https://img.youtube.com/vi/0gRVPLfZl7w/0.jpg)](https://www.youtube.com/watch?v=0gRVPLfZl7w)

* [Pretalx page](https://pretalx.com/juliacon2023/talk/RBHAER/)
* [Slides](https://www.docswell.com/s/hyrodium/5JL8EQ-JuliaCon2023)

## TL;DR
You can make a *holdable* smooth surface model with this repository.

![](docs/src/img/overview.png)

The main part of this project is how to determine a planer shape from a strip on the target curved surface.
In mathematics, this mapping is called "embedding".
We determined the embedding by minimizing its elastic strain energy.
This is the meaning of "Elastic Surface Embedding".

## Overview: How to make a surface model
### step 1: Define a shape of a surface (and split it into strips)
The definition must consist of parametric mapping and its domain.
For example, a paraboloid can be parametrized as below.

$$
\begin{aligned}
\boldsymbol{p}_{[0]}(u^1,u^2)
&= \begin{pmatrix}
u^1 \\
u^2 \\
(u^1)^2+(u^2)^2
\end{pmatrix} &
(u^1, u^2)
\in [-1,1] \times [-1,1]
\end{aligned}
$$

The domain will be split into $D^{(i)}$.

$$
\begin{aligned}
D^{(i)}
= [-1,1] \times \left[\frac{i-1}{10}, \frac{i}{10}\right] \qquad (i = 1,...,10)
\end{aligned}
$$

### step 2: Numerical analysis
This is the main part.
Split the surface into strips, and compute the embeddings.

```julia
using ElasticSurfaceEmbedding
using IntervalSets
using StaticArrays
# Overload the shape definition
ElasticSurfaceEmbedding.surface(x,y) = SVector(x, y, x^2+y^2)
# (1) split the surface into strips
dom = [(-1..1, (i-1)/10..i/10) for i in 1:10]
# (2) Embed the strips onto a plane
res = auto_allsteps(dom)
export_pinned_steps("paraboloid", res)
```

For more information, read [this document](https://hyrodium.github.io/ElasticSurfaceEmbedding.jl/dev/run-julia/).
The image below is a result for the domain $D^{(1)}$.

<img src="docs/src/img/bspline_strain.png" width="800">

### step 3: Edit on your favorite vector graphics editor
The output files are in SVG format.
After editing the SVG files, you can print the graphics or cut papers with a laser cutting machine.

<img src="docs/src/img/inkscape.png" width="800">

### step 4: Craft a paper model
This is the final step.
Cut papers into strips, and weave them into the surface.

<img src="docs/src/img/assembling.png" width="800">


## Directions: If you like..
### ..making crafts :scissors:
| <img src="docs/src/img/craft.png" align="top" height="150" width="150"> | Print Appendix B from [my paper on arXiv](https://arxiv.org/abs/2211.06372), and <a href="https://hyrodium.github.io/ElasticSurfaceEmbedding.jl/dev/craft/">make your own surface model. <br> Laser cutting machine is useful, but it's not necessary. |
| --- | :-- |

### ..computing :octocat:
| <img src="docs/src/img/juliawolfram.png" align="top" height="150" width="150"> | Clone this repository, and run the [Julia code](https://hyrodium.github.io/ElasticSurfaceEmbedding.jl/dev/run-julia/) or [Wolfram code](https://github.com/hyrodium/ElasticSurfaceEmbedding-wolfram)! <br> Any issues and pull requests are welcomed. |
| --- | :-- |

### ..mathematics or physics :globe_with_meridians:
| <img src="docs/src/img/math.png" align="top" height="150" width="150"> | Read [our paper on arXiv](https://arxiv.org/abs/2211.06372). Here's our theoretical framework: <br> ・Mathematical model: [Nonlinear elasticity](https://www.sciencedirect.com/topics/engineering/geometric-nonlinearity) on [Riemannian manifold](https://en.m.wikipedia.org/wiki/Riemannian_manifold) <br> ・Geometric representation: [B-spline manifold](https://en.wikipedia.org/wiki/Non-uniform_rational_B-spline) <br> ・Numerical analysis: [Galerkin method](https://en.wikipedia.org/wiki/Galerkin_method), [Newton-Raphson method](https://en.wikipedia.org/wiki/Newton%27s_method) |
| --- | :-- |

### ..me! :bowtie:
| <img src="docs/src/img/me.jpg" align="top" height="150" width="150"> | Follow [my twitter account](https://twitter.com/Hyrodium). <br> Visit [my website](https://hyrodium.github.io/). <br> Read [my paper on arXiv](https://arxiv.org/abs/2211.06372). |
| --- | :-- |

## Gallery
<img src="docs/src/img/Paraboloid1.png" width="160"> <img src="docs/src/img/Paraboloid2.png" width="160"> <img src="docs/src/img/Paraboloid3.png" width="160"> <img src="docs/src/img/Paraboloid4.jpg" width="160"> <img src="docs/src/img/Paraboloid5.png" width="160">

<img src="docs/src/img/CatenoidHelicoid.gif" width="400">

<img src="docs/src/img/stereographicprojection.png" width="800">
