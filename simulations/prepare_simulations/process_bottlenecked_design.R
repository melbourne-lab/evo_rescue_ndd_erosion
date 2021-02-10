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
  # Now merge with some population labels
  # alpha and evolve merge with sampno
  merge(y = data.frame(
    alpha  = rep(c(0, 0.0035), reps = 24*4),
    evolve = rep(c(TRUE, FALSE), each = 24*4)
  ) %>% 
    mutate(sampno = 1:192),
  by = 'sampno') %>%
  arrange(sampno, i)

# Test to make sure this worked correctly
pops.design1 %>%
  group_by(sampno) %>%
  mutate(n.pop0 = n()) %>%
  distinct(sampno, n.pop0, alpha, evolve) %>%
  group_by(n.pop0, alpha, evolve) %>%
  summarise(n.samps = n())

# Ah... there's a little bit of selection happening to get into certain design groups...  


