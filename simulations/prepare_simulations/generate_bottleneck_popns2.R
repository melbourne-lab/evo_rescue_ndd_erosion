# Adapted from script 'generate_bottlenecked_popns.R'
# Adapted and run by SN, July 6 2020

##### Load packages, namespace, etc. 

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

rm(list = ls())

setwd('~/demgen_allelic')

source('base_model_source/sim_functions.R')

# Initialize seed
set.seed(87795)

##### Parameters

n.trials = 5000
n.loci = 25
theta = 0
wfitn = sqrt(1 / 0.14 / 2)

##### Initialize populations

# This script will be more efficient if populations are initialized to all have
# one individual in each sex. That way, we don't waste trials on populations
# with two individuals of the same sex (who can not sustain a population).

# To do this, initialize one large population.
# Then, manually reset the sex column.
# This should still produce a genetically-representative population sample.
pop.large = init.sim(a = c(1/2, -1/2),
                     params = data.frame(init.row = 2 * n.trials,
                                         n.loci = n.loci, 
                                         # 2 individuals in each trial
                                         n.pop0 = 2 * n.trials,
                                         # Is used for calculating phenotypes and
                                         # number of offspring.
                                         w.max = 2, 
                                         theta = theta, 
                                         wfitn = wfitn,
                                         sig.e = 0.5,
                                         alpha = 0))

# Now: reset sex
# Resetting sex will also require resetting the number of offspring (r_i).
pop.large = pop.large %>%
  # Reset sex
  mutate(fem = rep(c(TRUE, FALSE), times = n.trials)) %>%
  # Recalculate fitness
  # (NOTE: setting alpha = 0 for speed)
  mutate(r_i = rpois(lambda = ifelse(fem, 2 * w_i, 0),
                     n = nrow(.)))

# Finally: set trial number
pop.large$trial = rep(1:n.trials, each = 2)

##### Run simulations

liszt = vector('list', n.trials)

# Note: the seed was set above when initializing populations.

for (trial.no in 1:n.trials) {
  liszt[[trial.no]] = sim(
    a = c(1/2, -1/2),
    params = data.frame(end.time = 9,
                        init.row = 1e4,
                        n.loci = n.loci, 
                        w.max = 2, 
                        theta = theta, 
                        wfitn = wfitn,
                        sig.e = 0.5,
                        alpha = 0),
    init.popn = pop.large %>%
      filter(trial %in% trial.no) %>%
      select(-trial),
    evolve = TRUE
  )
  print(trial.no)
}

all.trials = unroller(liszt)

##### Do the filtering (get populations reaching size 1000)

# Set threshold size above which to pick out populations.
thresh = 1000

# Pick out only the populations which reach size `thresh`
all.k.pops = all.trials %>%
  group_by(trial, gen) %>%
  mutate(n = n()) %>%
  group_by(trial) %>%
  filter(any(n > thresh)) %>%
  ungroup() %>%
  select(-n)

# Write this to a file, noting in name that all generations are here.
# (this could be useful for looking at loss of variance and the drift and
# fixation process!)
write.csv(all.k.pops, 'simulations/prepare_simulations/outputs/all_k_pops_all_gens.csv',
          row.names = FALSE)

# A more space-efficient data frame: has only the first generation at the given
# threshold size.
all.k.pops %>%
  group_by(trial, gen) %>%
  mutate(n = n()) %>% group_by(trial) %>%
  filter(gen == min(gen[n > thresh])) %>%
  ungroup() %>%
  select(-n) %>%
  write.csv(file = 'simulations/prepare_simulations/outputs/all_k_pops_min_gen.csv',
            row.names = FALSE)

# For fun: track the environmental variation over time.
var.k = all.k.pops %>%
  select(-c(g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(trial, i, gen)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(trial, gen, locus) %>%
  summarise(p = mean(val > 0),
            n = length(unique(i))) %>%
  group_by(trial, gen) %>%
  summarise(var.t = sum(2 * p * (1 - p) / 25),
            n.t = n[1])

write.csv(var.k, 'simulations/prepare_simulations/outputs/all_k_pops_var_loss.csv',
          row.names = FALSE)
