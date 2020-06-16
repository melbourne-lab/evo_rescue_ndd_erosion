# A script for generatin bottlenecked populations, run on melb2 for speed.
# Here, initialize population at size two, run a lot of these with the hopes
# that some will reach size 1000 at some point. These 1000-size populations will
# be the subsampled in simulations to get low-variance populations.

# NOTE: all of these thus far (june 16) have *no environmetal variance*.

##### Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

##### Clear namespace
rm(list = ls())

##### Set working directory
setwd('~/demall')

##### Load in source functions
source('base_model_source/sim_functions.R')

##### Parameters

# Initial trials suggested that populations initialized at size 2 take at least
# nine generations to reach a size of 1000. Additionally, running 500 trials
# yielded eight popuations which produced 1000 individuals by generation 9 (note
# that many of these populations were initialized with two males or two females
# and were incapable of naturally reproducing).

# I'll run 2000 trials here.

# Use 25 loci for now (initial trials suggest this is large enough)
# theta = 0 - populations here are already adapted to their population
# wfitn is the same value used in other situations

n.trials = 2000
n.loci = 25
theta = 0
wfitn = sqrt(1 / 0.14 / 2)

# Use threshold of 1000 (per proposal)
thresh = 1000

##### Run simulations

# Store sims in list object
liszt = vector('list', n.trials)

# Set seed
set.seed(87795)

# Run through trials
# Print trial number
for (trial in 1:n.trials) {
  liszt[[trial]] = sim(
    a = c(-1/2, 1/2),
    params = data.frame(end.time = 8,
                        init.row = 1e5,
                        n.loci = n.loci, 
                        n.pop0 = 2,
                        w.max = 2, 
                        theta = theta, 
                        wfitn = wfitn,
                        sig.e = 0)
  )
  print(trial)
}

# Put trials into one (very large!) data frame.
all.trials = unroller(liszt)

# For space/memory, pick out only the trials with one parent of each sex.
init.pops = all.trials %>%
  group_by(trial) %>%
  filter(mean(fem[gen %in% 1]) == 0.5)

##### Do the filtering

# Pick out only the populations which reach size `thresh`
all.k.pops = init.pops %>%
  group_by(trial, gen) %>%
  mutate(n = n()) %>%
  group_by(trial) %>%
  filter(any(n > thresh)) %>%
  ungroup() %>%
  select(-n)

# Write this to a file, noting in name that all generations are here.
# (this could be useful for looking at loss of variance and the drift and
# fixation process!)
write.csv(all.k.pops, 'bottlenecks/out/all_k_pops_all_gens.csv',
          row.names = FALSE)


# A more space-efficient data frame: has only the first generation at the given
# threshold size.
all.k.pops %>%
  group_by(trial, gen) %>%
  mutate(n = n()) %>% group_by(trial) %>%
  filter(gen == min(gen[n > thresh])) %>%
  ungroup() %>%
  select(-n) %>%
  write.csv(file = 'bottlenecks/out/all_k_pops_min_gen.csv',
            row.names = FALSE)

