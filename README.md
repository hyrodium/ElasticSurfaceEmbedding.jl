# Elastic Surface Embedding; Weaving Parer Strips

## TL;DR
You can make a *holdable* smooth surface model with this repository.

![](docs/src/img/overview.png)

The main part of this project is how to determine a planer shape from a strip on curved surface.
In mathematics, this mapping is called "embedding".
We determined the embedding by minimizing its elastic strain energy.
This is the meaning of "Elastic Surface Embedding".

## Overview: How to make a surface model
### step 1 : Define a shape of surface (and split into strips)
The definition must consists of parametric mapping and its domain.
For example, a paraboloid can be parametrized as below.

<img src="docs/src/img/paraboloid-formula.png" width="325">

The domain D will be split into D_i.

<img src="docs/src/img/paraboloid-domain.png" width="465">

### step 2 : Numerical analysis
This is the main part.
Split the surface into pieces, and compute the Eucledian embedding.
For more information, read [this document](/Julia).
The image below is a result for the domain D_1.

<img src="docs/src/img/NurbsStrain.png" width="800">

### step 3 : Edit on vector graphics editor
The output files are SVG format.
After editing the svg files, you can print the graphics or cut papers by laser cutting machine.

<img src="docs/src/img/inkscape.png" width="800">

### step 4 : Craft a paper model
This is the final step.
Cut papers into strips, and weave them into surface.

<img src="docs/src/img/assembling.png" width="800">


## Directions: If you like..
### ..making crafts:scissors:
| <img src="docs/src/img/craft.png" align="top" height="150"> | Download the [Paraboloid example](/Example/Paraboloid.pdf) and [make your own surface model](Craft). <br> Laser cutting machine is useful, but it's not necessary. |
| --- | :-- |

### ..computing:octocat:
| <img src="docs/src/img/juliawolfram.png" align="top" height="150"> | Clone this repository, and run the [Julia code](/Julia) or [Wolfram code](/Wolfram)! <br> Any issues and pull requests are welcomed. |
| --- | :-- |

### ..mathematics or physics:globe_with_meridians:
| <img src="docs/src/img/math.png" align="top" height="150"> | Read our upcoming paper. Here's our theoretical framework: <br> ・Mathematical model: [Nonlinear elasticity](https://www.sciencedirect.com/topics/engineering/geometric-nonlinearity) on [Riemannian manifold](https://en.m.wikipedia.org/wiki/Riemannian_manifold) <br> ・Geometric representation: [NURBS (esp. B-spline manifold)](https://en.wikipedia.org/wiki/Non-uniform_rational_B-spline) <br> ・Numerical analysis: [Galerkin method](https://en.wikipedia.org/wiki/Galerkin_method), [Newton-Raphson method](https://en.wikipedia.org/wiki/Newton%27s_method) |
| --- | :-- |

### ..me!:bowtie:
| <img src="docs/src/img/me.jpg" align="top" height="150"> | Follow [my twitter account](https://twitter.com/Hyrodium). <br> Visit [my website](https://hyrodium.github.io/). <br> Read my upcoming paper. |
| --- | :-- |

## Gallery
<img src="docs/src/img/Paraboloid1.png" width="160"> <img src="docs/src/img/Paraboloid2.png" width="160"> <img src="docs/src/img/Paraboloid3.png" width="160"> <img src="docs/src/img/Paraboloid4.jpg" width="160"> <img src="docs/src/img/Paraboloid5.png" width="160">

<img src="docs/src/img/CatenoidHelicoid.gif" width="400">
