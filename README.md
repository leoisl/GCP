A header-only c++ implementation of genotype confidence percentile (GCP) production.

## Usage

Copy "GCP.h" in your project.

The requirements are:

1. Defining a struct or class holding data: let's call it `ModelData`.
1. Defining a genotyper class. This **must have**:
	* A constructor taking a `ModelData` object.
	* A method called `get_genotype_confidence`, with signature `double get_genotype_confidence()`,
		returning the confidence of the genotype call.
1. Overriding a method to simulate data produced as a `ModelData`.

Only point 3 is actually required by GCP. Points 1 and 2 would usually be internally implemented as part of a tool's genotyping.

Follow this [example](examples/Full.cpp) for a concrete walkthrough!

## Users

* [pandora][pandora]
* [gramtools][gramtools]

## Credit

Originally based on minos Genotype confidence percentile (https://github.com/iqbal-lab-org/minos/blob/ed975d1099c6403eff79b04b0d2064eebfa95e73/minos/genotype_confidence_simulator.py) by Martin Hunt


[gramtools]: https://github.com/iqbal-lab-org/gramtools
[pandora]: https://github.com/rmcolq/pandora
