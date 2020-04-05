# init_model

#### What is here

This folder contains the original model code and some examples from March, 2020 model development.

This model features many (but not all) of the same demographic features from the original demographic and genetic model, meant to model evolutionary rescue. However, the genetic details are different. Rather than following the genetic framework as outlined by Gomulkiewicz and Holt, 1995 and Lande, 1976, this model features explicit alleles which determine genotypic values. This allows a mechanistic and dynamic genetic variance.

Because the original model was meant solely to look at change in genetic diversity over time, this model is missing a couple of features and/or has some discrepancies.

- Fitness; the Gaussian fitness term here is exp(- (z - theta)^2 / wsq); in Lande and Gomulkiewcz and Holt, this denominator is (2 * w^2); simulation trials were run with (1 / 0.14) for wsqd which is not the same value used in previous simulations
- Density dependence; no Ricker term in this model
- Non-heritable contributions to phenotype; this wasn't incorporated yet (i.e., phenotype is entirely determined by genotype)

#### Directory contents

*sim_evolve_allele_spec_brs.R* - holds main simulation functions and some auxiliary analysis functions

*test_allele_spec_brs.R* - runs 1000 simulations for one allele frequency; exports outputs

*example_evolve_allele_spec.R* - useful example script for running through one simulation instance, plus code, and also running through summaries of the 1000 sims run above

*server_calc_sima_a.R* - script run on server for calculating additive genetic variance in 1000k trials

*analyze_geno_change.R* - script for analyzing genotypic change in 1000 simulation runs (crudely)

*eg_1k_n0_20.csv* - summary statistics for 1000 simulation runs

*allelic_all_sigma_a.csv* - summary genetic variance for 1000 simulation runs

*notes_2020-03-18.md* - meeting notes from initial analysis of model

