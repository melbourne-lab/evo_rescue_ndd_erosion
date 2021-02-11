# Generate bottlenecks for testing experimental designs of Triboleum experiments.
# Run on server - SN 9 Feb 2021

# Load packages
library(tidyselect)
library(dplyr)
library(tidyr)

# Clear namespace
rm(list = ls())

# Load functions for sims and processing
source('base_model_source/sim_functions.R')

# Number of trials to prepare for (n.b. fewer run)
n.trials = 200

# Parameters
# n.b., 2 individuals for each simulation
# n.b., w.max influences the rate of population growth, but also, the
# spread of individual phenotypes in the final populations
# n.b., theta = 0 populations are adapted to their environments
pars = data.frame(
  init.row = 2e4,
  end.time = 10,
  n.pop0 = 2 * n.trials,
  n.loci = 25,
  w.max = 2,
  theta = 0,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0
)

# Set seed
set.seed(1004)

# Take the initialized population - 
# reassign sex to be female/male (prevents immediate extincton)
# try - trial (two individuals per trial only)
# reassign fecundities (to correspond to reassigned sex)
sample.from = init.sim(params = pars) %>%
  mutate(fem = rep(c(T, F), times = n.trials),
         try = rep(1:n.trials, each = 2),
         i   = rep(1:2, times = n.trials)) %>%
  mutate(r_i = rpois(n = 2*n.trials, lambda = ifelse(fem, 2 * w_i, 0)))

# head(sample.from)

liszt = vector('list', n.trials)

# Simulate 220 times (note: population is initialized to two parents undr `try`)
for (trial in 1:n.trials) {
  
  neck1 = propagate.sim(popn = sample.from %>% 
                          filter(try %in% trial) %>% 
                          select(-try), 
                        params = pars)
  
  if (nrow(neck1) > 1) {
    neck2 = neck1 %>% 
      sample_n(2) %>%
      mutate(fem = c(TRUE, FALSE)) %>%
      mutate(r_i = rpois(n = 2, lambda = ifelse(fem, 2 * w_i, 0)))
    
    if (sum(neck2$r_i)) {
      liszt[[trial]] = sim(params = pars, init.popn = neck2)
      print(trial)
    } else {
      print('fail step 2')
    }
    
  } else {
    print('fail step 1')
  }
  
}

# Save the trials which fail (can go back and compare to original populations)

survivors = liszt[sapply(liszt, function(x) !is.null(x))]

sapply(survivors, function(x) x %>% filter(gen %in% max(gen)) %>% nrow())

# Subsetting this liszt out in a savvy way
# (because we will be dealing an inordinate amount of junk data otherwise)

# Wanted: 24 source populations, each with 
# (20 + 100) * 2 * 2 individuals
# (2 population sizes, 2 evolution treatments, 2 alpha vals.)

# Here: give me indices of first population which are sufficent size for design 1
sum(sapply(survivors, function(x) x %>% filter(gen %in% max(gen)) %>% nrow()) > 480)
# At current seed value we can do it here.

surv.subset = which(
  sapply(survivors, function(x) x %>% filter(gen %in% max(gen)) %>% nrow()) > 480
)[1:24]

all.data = unroller(survivors[surv.subset])

# Write a file (with only the last generation in each sim)
all.data %>%
  group_by(trial) %>%
  filter(gen %in% max(gen)) %>%
  ungroup() %>%
  write.csv(file = 'experiment/expmt_bottleneck24/out/source_pops_large.csv',
            row.names = FALSE)
