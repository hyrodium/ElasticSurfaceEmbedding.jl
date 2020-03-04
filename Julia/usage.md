# How to use

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
    * JLD
    * ForwardDiff
    * DifferentialEquations
    * FastGaussQuadrature
    * JSON
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
[flow chart here]

1. define the shape
1. initial



