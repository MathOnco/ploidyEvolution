# S-phase duration as a function of dNTP substrate concentration

Using ploidy as proxy of the number of Polymerases expressed in a cell and prior estimates of Michaelis-Menten parameters for the average DNA replication rates of a single polymerase, we project how the time required for a cell to complete S-phase increases with decreasing dNTP concentrations in stomach cancer (STAD) and with decreasing dNTP substarte concentration in Glioblastoma (GBM).

### Prerequisites

To run the code you will need to install:

```
R (version 3.6 or greater)
R packages:
        RMySQL
	xlsx
	cloneid
	devtools
	FactoMineR
	factoextra
	scatterplot3d
	akima
	flexclust
	matlab
	fields
	rgl
	plyr
```

### Installing

```
1. Start R and add the required packages
2. Load R function via "source('EnergyParameterization_dNTP_O2.R')" 
```

## Running
To run the function, type

```
EnergyParameterization_dNTP_O2(cancer = "GBM")
EnergyParameterization_dNTP_O2(cancer = "STAD")
```


