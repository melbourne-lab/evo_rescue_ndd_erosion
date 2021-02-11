# Process bottlenecked populations for simulating our experimental design
# SN - February 10 2021

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Clear namespace
rm(list = ls())

# Read in final gen in simulated populations
simpops = read.csv('simulations/prepare_simulations/outputs/bottlenecks_for_design.csv')

rawpops = simpops

simpops = rawpops

# Ensure this was done correctly
head(simpops)

# How large is each population?
table(simpops$trial)

### Preparing design 1
# Design 1 is one bottlenecked population, subsampled enough to do all trials
# Total number of trials is 24 * 4 (24 reps * 2 evolving treatments * 2 alpha trts)
# at each population size (20, 100).
# So, we want a population large enough to distinctly sample 96*100 and 96*20
# Pick out these indices.
inds.design1 = which(table(simpops$trial) > (96*100 + 96*20))[1]

pops.design1 = simpops %>% filter(trial %in% inds.design1)

pops.design1 = pops.design1 %>%
  # Get only as many individuals as we need:
  slice(1:(96*100 + 96*20)) %>%
  # Designate design number
  mutate(design = 1) %>%
  # Assign some labels for population sizes
  # `sampno` is a sample number (to differentiate samples)
  mutate(i = rep(c(1:20, 1:100), times = 96),
         idiff = c(0, diff(i)),
         sampno = 1 + cumsum(idiff < 0)) %>%
  select(-idiff) %>%
  # Add population size variable
  group_by(sampno) %>%
  mutate(n.pop0 = max(i)) %>%
  # Now merge with some population labels
  # alpha and evolve merge with sampno
  # Note: due to indexing above, n.pop0 must be first!
  merge(y = expand.grid( 
    n.pop0 = c(20, 100),
    alpha  = c(0, 0.0035),
    evolve = c(TRUE, FALSE),
    block  = 1:4,
    rep    = 1:6
  ) %>% 
    mutate(sampno = 1:192), 
  by = c('sampno', 'n.pop0'))

# Test to make sure this worked correctly
pops.design1 %>%
  distinct(sampno, n.pop0, alpha, evolve, block, rep) %>%
  group_by(n.pop0, alpha, evolve, block, rep) %>%
  summarise(n.samps = n()) %>%
  spread(key = rep, value = n.samps)

pops.design1 %>%
  distinct(sampno, n.pop0, alpha, evolve, block, rep) %>%
  group_by(n.pop0, alpha, evolve, block, rep) %>%
  summarise(n.samps = n()) %>%
  ungroup() %>%
  summarise(npops = sum(n.samps))

# Ah... there's a little bit of selection happening to get into certain design groups...  

### Preparing design 2
# Design 2 is one bottlenecked population, recycled/subsampled for each block
# Total number of trials is 6 * 4 (4 reps * 2 evolving treatments * 2 alpha trts)
# at each population size (20, 100).
# So, we want a population large enough to distinctly sample 96*100 and 96*20
# Pick out these indices.
inds.design2 = which(table(simpops$trial) > (24*100 + 24*20))[1:4]

pops.design2 = simpops %>% filter(trial %in% inds.design2)

pops.design2 = pops.design2 %>%
  group_by(trial) %>%
  # Get only as many individuals (per trial) as we need:
  slice(1:(24*100 + 24*20)) %>%
  # Ungroup to get sampno up to 192
  ungroup() %>%
  # Designate design number
  mutate(design = 2) %>%
  # Assign some labels for population sizes
  # `sampno` is a sample number (to differentiate samples)
  mutate(i = rep(c(1:20, 1:100), times = 96),
         idiff = c(0, diff(i)),
         sampno = 1 + cumsum(idiff < 0)) %>%
  select(-idiff)  %>%
  # Add population size variable
  group_by(sampno) %>%
  mutate(n.pop0 = max(i)) %>%
  # Now merge with some population labels
  # alpha and evolve merge with sampno
  # Note: due to indexing above, n.pop0 must be first, block needs to be last!!!
  merge(y = expand.grid(
    n.pop0 = c(20, 100),
    alpha  = c(0, 0.0035),
    evolve = c(TRUE, FALSE),
    rep    = 1:6,
    block  = inds.design2
  ) %>% 
    mutate(sampno = 1:192),
  by.x = c('trial', 'sampno', 'n.pop0'),
  by.y = c('block', 'sampno', 'n.pop0')) %>%
  arrange(sampno, i)
  
# Test to make sure this worked correctly
pops.design2 %>%
  distinct(sampno, n.pop0, alpha, evolve, trial, rep) %>%
  group_by(n.pop0, alpha, evolve, trial, rep) %>%
  summarise(n.samps = n()) %>%
  spread(key = rep, value = n.samps)

pops.design2 %>%
  distinct(sampno, n.pop0, alpha, evolve, trial, rep) %>%
  group_by(n.pop0, alpha, evolve, trial, rep) %>%
  summarise(n.samps = n()) %>%
  ungroup() %>%
  summarise(npops = sum(n.samps))
  
### Preparing design 3
# Design 3 is one bottlenecked population per unique trial
# Total number of trials is 96
# (6 reps * 4 blocks * 2 evolving treatments * 2 alpha trts * 2 size treatments)
# Do I have enough trials for that?
length(unique(simpops$trial))
# Oh shoot I messed up and only did 96... probably greater than 120.
# Well I guess this will be a paired design, with pairs of large and small.
all(table(simpops$trial) > 120)
# Use all 96 trials (i.e., no need to subset simpops)

pops.design3 = simpops %>%
  group_by(trial) %>%
  # Get only as many individuals (per trial) as we need:
  slice(1:120) %>%
  # Ungroup to get sampno up to 192
  ungroup() %>%
  # Designate design number
  mutate(design = 3) %>%
  # Assign some labels for population sizes
  # `sampno` is a sample number (to differentiate samples)
  mutate(i = rep(c(1:20, 1:100), times = 96),
         idiff = c(0, diff(i)),
         sampno = 1 + cumsum(idiff < 0)) %>%
  select(-idiff)  %>%
  # Add population size variable
  group_by(sampno) %>%
  mutate(n.pop0 = max(i)) %>%
  # Now merge with some population labels
  # alpha and evolve merge with sampno
  # Note: due to indexing above, n.pop0 must be first
  merge(y = expand.grid(
    n.pop0 = c(20, 100),
    alpha  = c(0, 0.0035),
    evolve = c(TRUE, FALSE),
    rep    = 1:6,
    block  = 1:4
  ) %>% 
    mutate(sampno = 1:192),
  by = c('sampno', 'n.pop0')) %>%
  arrange(sampno, i)

pops.design3 %>%
  distinct(sampno, n.pop0, alpha, evolve, block, rep) %>%
  group_by(n.pop0, alpha, evolve, block, rep) %>%
  summarise(n.samps = n()) %>%
  spread(key = rep, value = n.samps)

pops.design3 %>%
  distinct(sampno, n.pop0, alpha, evolve, block, rep) %>%
  group_by(n.pop0, alpha, evolve, block, rep) %>%
  summarise(n.samps = n()) %>%
  ungroup() %>%
  summarise(npops = sum(n.samps))

nrow(pops.design3)

### Write csvs

# Design 1
write.csv(pops.design1, row.names = FALSE,
          'simulations/prepare_simulations/outputs/bottlenecks_design1.csv')

# Design 2 (going to rename one column)
write.csv(pops.design2 %>% rename(block = trial), row.names = FALSE,
          'simulations/prepare_simulations/outputs/bottlenecks_design2.csv')

# Design 3
write.csv(pops.design3, row.names = FALSE,
          'simulations/prepare_simulations/outputs/bottlenecks_design3.csv')

