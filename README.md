# Allele-specific models of evolutionary rescue

Simulation and analysis code for the manuscript "How density dependence, genetic erosion and the extinction vortex impact evolutionary rescue" by Nordstrom et al. 

Requirements:
- R (version used: 4.1.2)
- dplyr (1.0.7)
- tidyr (1.1.13)
- rstanarm (2.21.1)

For figures:
- ggplot2 (3.3.4)
- cowplot (1.1.1)

## Repository structure

`sim_functions.R`: Source script for simulation functions and other handy wrappers

Directory `simulations`:

* `simulations/run_final_scripts/` scripts for running simulations used in analysis
	
	+ `simulations/run_final_scripts/sim_alldata_n{20/100}_a{000/035}_{hi/low}var.R`: script for running 15-generation simulations with size 20/100, alpha 0/0.0035, high/low genetic diversity
		
	+ `simulations/run_final_scripts/longsims_n{20/100}_a{000/035}_{hi/low}var.R`: script for running 50-generation simulations with size 20/100, alpha 0/0.0035, high/low genetic diversity
		
	+ `simulations/combine_alldata.R`: short script for concatenating 15-generation simulation output files into one output file
		
* `simulations/outputs/` subdir for holding CSV outputs from simulations (used in analysis)
	
	+ `simulations/outputs/alldata_n{20/100}_a{000/035}_{hi/low}var.csv`: 15-generation simulation output (size 20/100, alpha 0/0.0035, high/low genetic diversity)
		
	+ `simulations/outputs/longsims_n{20/100}_a{000/035}_{hi/low}var.R`: 50-generation simulation output (size 20/100, alpha 0/0.0035, high/low genetic diversity)

Directory `analysis_results`:

* Contains individual scripts labelled to reflect contents (including tables, figures, and supplemental information)

## Contact information

Scott Nordstrom (scottwatsonnordstrom&gmail)

## Changelog

Readme finalized October 24, 2022.
