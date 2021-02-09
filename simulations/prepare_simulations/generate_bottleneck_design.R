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
n.trials = 400

# Parameters
# n.b., 2 individuals for each simulation
# n.b., w.max influences the rate of population growth, but also, the
# spread of individual phenotypes in the final populations
# n.b., theta = 0 populations are adapted to their environments
pars = data.frame(
  init.row = 2e4,
  end.time = 15,
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
         try = rep(1:n.trials, each = 2)) %>%
  mutate(r_i = rpois(n = 2*n.trials, lambda = ifelse(fem, 2 * w_i, 0)))

# head(sample.from)

# Initialize 220 trials
# (an inital run of this script suggested 201 simulations is sufficient)
liszt = vector('list', 220)

# Simulate 220 times (note: population is initialized to two parents undr `try`)
for (trial in 1:220) {
  liszt[[trial]] = sim(params = pars, init.popn = sample.from %>% filter(try %in% trial))
  print(trial)
}

# Subsetting this liszt out in a savvy way
# (because we will be dealing an inordinate amount of junk data otherwise)

# Threshold 1: one source population
# needs to be met once
# actually you know just slice and dice four times (this is easier and will be quicker)
thresh1 = sum(rep(c(20, 100), each = 24*4)) # 24 reps for each init. size x alpha x evolving trt

# Threshold 2: four source populations (one per block)
# needs to be met six times
thresh2 = sum(rep(c(20, 100), each = 4*4))
times2 = 6

# Threshold 3: 24 source populations (one per trial)
# eh... meet this 4*24 = 96  times (once per combo)
thresh3 = 120 # 100 + 20
times3 = 96

# Print for check
print('thresh assigned')

# Here: give me indices of first population which are sufficent size for design 1
set1 = which(
  sapply(liszt, function(x) x %>% filter(gen %in% max(gen)) %>% nrow()) > thresh1
)[1]

# Here: give me indices of first (times2) populations which are sufficent size
# for design 2
set2 = which(
  sapply(liszt, function(x) x %>% filter(gen %in% max(gen)) %>% nrow()) > thresh2
)[1:times2]

# Here: give me the indices of the first (times3) populations which are
# sufficient size for design 3.
set3 = which(
  sapply(liszt, function(x) x %>% filter(gen %in% max(gen)) %>% nrow()) > thresh3
)[1:times3]

print('sets designated')

print(set1)
print(set2)
print(set3)

# Subset out only the trials indicated above, cat them together
all.data = unroller(liszt[unique(c(set1, set2, set3))])

# Write a file (with only the last generation in each sim)
all.data %>%
  group_by(trial) %>%
  filter(gen %in% max(gen)) %>%
  ungroup() %>%
  write.csv(file = 'simulations/prepare_simulations/outputs/bottlenecks_for_design.csv',
            row.names = FALSE)
