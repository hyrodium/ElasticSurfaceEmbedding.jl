# Julia code for Elastic Surface Embedding

## Overview
### Input
* parametric mapping
* its domain

### Output
* Slack export
* (Result directory)
  * Slack
  * strain
  * ..

## Prerequisite
In this repository, Julia is for numerical analysis of surface embedding.

* [Julia](https://julialang.org/) (≥v1.0)
* Packages
    * IntervalSets
    * Luxor
    * Colors
    * ForwardDiff
    * DifferentialEquations
    * FastGaussQuadrature
    * JSON
	* Images
	* ImageMagick (for linux)
* Slack (optional)
    * Incoming webhook URL
    * OAuth access token
    * Channel ID

## JSON settings
Before running code, you should configure `config.json`.
Here's an example:

```
{
	"Poisson's Ratio": 0.25,
	"Default Number of Integration Points": 10,
	"Output Directory": "~/ElasticSurfaceEmbedding-Result/",
	"Slack":{
		"Incoming Webhook Url": "https://hooks.slack.com/services/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",
		"OAuth Access Token": "xoxp-587883028310-587883028694-622674491761-a37533ee7606e19d2216a6a88f185641",
		"Channel ID": "CJK5ALKJS"
	}
}
```

## Running code
You can run the code in REPL, Juno, or Jupyter notebook.
I strongly reccomend Juno.

### Flowchart
The computation process proceeds like this flowchart:

[flow chart here]

* Define the shape of surface
* Split the surface into strips
* Check the strain prediction
* Initial configuration
* Newton-Raphson method itteration
* Refinement B-spline manifold
* Finish Computation

### (Load packages)
```
using IntervalSets
...
```

### Define the shape of surface
Through this document, we treat a paraboloid as an example.

```
@ParametricMapping ...
```

### Split the surface into strips
Define a domain of the strip with symbol `D`.

```
D = [0..1,0..1]
```

### Check the strain prediction

```
ShowStrain(D)
```

Negative number means compression, and positive number means tension.
The absolute value of these numbers should be less than 0.01.

### Initial configuration
If you finished checking the strain prediction, the next step is determine the initial configuration.

```
InitialConfiguration(D)
```


### Newton-Raphson method itteration

```
NewtonMethod(D)
```

You can choose the fixing method:
* hoge (default)
* fuga
* FixThreePoints


### Refinement B-spline manifold


```
Refinement()
```

### Finish

```
FinishComputation()
```