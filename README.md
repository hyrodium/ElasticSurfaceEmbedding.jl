# Elastic Surface Embedding

## Overview
You can make a smooth surface model by this code.

<img src="img/Paraboloid1.png" width="160"> <img src="img/Paraboloid2.png" width="160"> <img src="img/Paraboloid3.png" width="160"> <img src="img/Paraboloid4.jpg" width="160"> <img src="img/Paraboloid5.png" width="160">

## Prerequisite
### Julia
In this repository, Julia is for numerical analysis of surface embedding.

* [Julia](https://julialang.org/)(≥v1.3)
* Packages
    * IntervalSets
    * Luxor
    * Colors
    * JLD
    * ForwardDiff
    * DifferentialEquations
    * FastGaussQuadrature
    * JSON
* Slack (optional)
    * Incoming webhook URL
    * OAuth access token
    * Channel ID

### Wolfram Engine (optional)
The wolfram code is for proof of the main theorem of my theory.
So, you don't need part if you just want to make a surface model.

* [Wolfram Engine](https://www.wolfram.com/engine/)


## How to make a surface model
### Step 1 : Define shape


### Step 2 : Numerical analysis
This is the main part.
Split the surface into pieces, and compute the Eucledian embedding.
For more information, read [this](/Julia/usage.md).

### Step 3 : Edit on vector graphics editor
The output files are svg format.
You can print the graphics or cut papers by laser cutting machine, and this step is for preparation.

### Step 4 : Craft a paper model



## Theoretical framework

* Elasticity
    * Nonlinear elasticity (geometric non-linearity) on Riemannian manifold
* Numerical analysis
    * [NURBS](https://en.wikipedia.org/wiki/Non-uniform_rational_B-spline)
    * [Galerkin method](https://en.wikipedia.org/wiki/Galerkin_method)
    * [Newton-Raphson method](https://en.wikipedia.org/wiki/Newton%27s_method)

## More information
Visit [my website](https://hyrodium.github.io/Profile/) or read my upcoming paper.
