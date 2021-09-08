# Computing

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
* Slack (optional)
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
		"OAuth Access Token": "xoxb-1018669471475-1056240006178-dMLt4osKWMZ8QbQCirODnqtT",
		"Channel ID": "C010HSTF9LJ"
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
* Newton-Raphson method iteration
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
@parametric_mapping ...
```

### Split the surface into strips
Define a domain of the strip with symbol `D`.

```
D(i,n) = (-1.0..1.0, (i-1)/n..i/n)
```

### Check the strain prediction

```
i = 1
print_strain(D(i,10))
```

Negative number means compression, and positive number means tension.
The absolute value of these numbers should be less than 0.01.

### Initial configuration
If you finished checking the strain prediction, the next step is determine the initial configuration.

```
print_strain(D(i,10))
```


### Newton-Raphson method iteration

```
newton_onestep(fixingmethod=:FixThreePoints)
```

You can choose the fixing method from below:
* `:DefaultOrientation` (default)
* `:FixThreePoints`

### Refinement B-spline manifold


```
spline_refinement(p₊=[0,1],k₊=[Knots([]),Knots([(i-1/2)/10])])
```

### Pin the state

```
pin_state(tag="paraboloid-"*string(i+1))
```

### Finish

```
FinishComputation()
```

### Utilities


```
RemovePin(16)
```
```
computed_shapes()
```
```
print_knots()
```