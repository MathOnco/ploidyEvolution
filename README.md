# In-silico evolution of aneuploidy through chromosome mis-segregations

Mathematical models and pipelines used to infer model parameters from multi-modal data analysis. Sequencing data is analyzed to identify co-existing clones in a tumor sample; quantify their ploidy, copy number state and missegregations.

## Getting Started

Before being able to execute the code you will need to install Julia (≥1.5) and all the required dependencies.

### Prerequisites

To run the code you will need to install:

```
Julia (version 1.5 or greater)
Julia packages:
	KrylovKit
	LinearAlgebra
	DelimitedFiles
	DifferentialEquations
```

### Installing
On Mac:

```
1. Download Julia from https://julialang.org/downloads/
2. Open julia in terminal and add the required packages
3. Add package by typing ```]``` then ```add <package name>```
```

## Running tests
To run the script, open Julia from terminal and type

```
include("path/to/file/<filename>")
```

where ```filename``` is either ```computeCriticalCurveBisection.jl``` or ```testPloidyMovement_v4.jl```.

## Authors

* **Gregory Kimmel**
* **Philipp Altrock**
* **Noemi Andor**

## License

## Acknowledgements