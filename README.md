# Evolutionary rescue 

Simulation and analysis code for the paper "How density dependence, genetic erosion and the extinction vortex impact evolutionary rescue" by Nordstrom et al. 

Requirements:
- R (version used: 4.1.2)
- dplyr (1.0.7)
- tidyr (1.1.13)
- rstanarm (2.21.1)

For figures:
- ggplot2 (3.3.4)
- cowplot (1.1.1)

The robustness checks also use the parallel package (4.1.2).

## Repository structure

`sim_functions.R`: Source script for simulation functions and other handy wrappers
`sim_bh_functions.R`: Source script for simulations with Beverton-Holt density dependence (for robustness analysis)

Directory `simulations`:

* `simulations/run_final_scripts/` scripts for running simulations used in analysis
	
	+ `simulations/run_final_scripts/sim_alldata_n{20/100}_a{000/035}_{hi/low}var.R`: script for running 15-generation simulations with size 20/100, alpha 0/0.0035, high/low genetic diversity
		
	+ `simulations/run_final_scripts/longsims_n{20/100}_a{000/035}_{hi/low}var.R`: script for running 50-generation simulations with size 20/100, alpha 0/0.0035, high/low genetic diversity
		
	+ `simulations/combine_alldata.R`: short script for concatenating 15-generation simulation output files into one output file

* `simulations/run_robust_checks` scripts for running robust checks (one per parameter for each of `nloci` (number of loci - m in manuscript), `sig2e` (environmental contribution to phenotypic variance, sig^2_e in manuscript), `wfitn` (width of selection surface), `w_max` (maximum intrinsic fitness), `theta` (degree of environmental change), and `bholt` (with Beverton-Holt density dependence))
		
* `simulations/outputs/` subdir for holding CSV outputs from simulations (used in analysis)
	
	+ `simulations/outputs/alldata_n{20/100}_a{000/035}_{hi/low}var.csv`: 15-generation simulation output (size 20/100, alpha 0/0.0035, high/low genetic diversity)
		
	+ `simulations/outputs/longsims_n{20/100}_a{000/035}_{hi/low}var.R`: 50-generation simulation output (size 20/100, alpha 0/0.0035, high/low genetic diversity)

	+ subdirectory `robust_checks_out` with one output csv for each robustness check

Directory `analysis_results`:

* Contains individual scripts labelled to reflect contents (including tables, figures, and supplemental information)

* Contains subdirectory `robustness_checks` with code for making robustness check figures

## Contact information

Scott Nordstrom (scottwatsonnordstrom&gmail)

## Changelog

Readme finalized October 24, 2022.

Additional analysis for new submission finalized and merged in June 2, 2023.

Revisions for Figure 5 (legend addition) and source code for Beverton-Holt simulations added September 22, 2023.

Code for robustness checks (simulations and figure production) added October 13, 2023.

Readme updated October 14, 2023.

Correct and up-to-date Fig. 2 script (including analytical solutions) added to repo after being mistakenly left out, October 20, 2023.
