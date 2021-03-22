# In-silico evolution of aneuploidy through chromosome mis-segregations

Mathematical models and pipelines used to infer model parameters from multi-modal data analysis. Sequencing data is analyzed to identify co-existing clones in a tumor sample; quantify their ploidy, copy number state and missegregations.

## Getting Started

Before being able to execute the code you will need to install Julia (≥1.5) and all the required dependencies.

### Prerequisites

To run the code you will need to install:

```
Julia (version 1.5 or greater)
```

### Installing
On Mac, the following set of instructions will get all packages and corresponding dependencies needed to run the project.

```
1. Open terminal and go to directory where you want to clone repository
2. clone git https://github.com/MathOnco/ploidyEvolution.git
3. cd ploidyEvolution/Julia
4. julia --project=. install.jl
```

## Running tests
Run the script from terminal by

```
julia --project=. -L ploidyMovement.jl driver.jl
```

for help on running and options type

```
julia driver.jl -h
```


## Authors

* **Gregory Kimmel**
* **Philipp Altrock**
* **Noemi Andor**

## License

## Acknowledgements