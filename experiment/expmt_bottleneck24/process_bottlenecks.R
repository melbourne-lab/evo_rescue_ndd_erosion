# Process bottlenecked populations for simulating our experimental design
# SN - February 11 2021

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Clear namespace
rm(list = ls())

# Read in final gen in simulated populations
simpops = read.csv('experiment/expmt_bottleneck24/out/source_pops_large.csv')

rawpops = simpops

simpops = rawpops

# How large is each population? Should all be greater than 480.
table(simpops$trial)

# Now... do subsetting.

sourcepops = simpops %>%
  # For each trial, pick only the first 480 individuals
  group_by(trial) %>%
  slice(1:480) %>%
  # Now, relabel individuals
  mutate(i = rep(c(1:20, 1:100), times = 4),
         sampno = cumsum(i %in% 1)) %>%
  # Get population size per sample
  group_by(trial, sampno) %>%
  mutate(n.pop0 = n()) %>%
  ungroup() %>%
  # Merge in rest of experimental design
  # (assigning alpha, extinct to each sample)
  merge(y = expand.grid(
    n.pop0 = c(20, 100),
    evolve = c(TRUE, FALSE),
    alpha  = c(0, 0.0035)
  ) %>%
    mutate(sampno = 1:8),
  by = c('sampno', 'n.pop0')
  ) %>%
  arrange(trial, sampno, i)

nrow(sourcepops)

# Test to make sure this worked correctly
sourcepops %>%
  distinct(sampno, n.pop0, alpha, evolve, trial) %>%
  group_by(n.pop0, alpha, evolve) %>%
  summarise(n.samps = n())

# Rename `trial` to `sourcepop` and break that down into
# block and rep
sourcepops = sourcepops %>%
  rename(sourcepops = trial) %>%
  mutate(block = (sourcepops-1) %%  4 + 1,
         repct = (sourcepops-1) %/% 4 + 1)

# Export
write.csv(sourcepops, row.names = FALSE,
          file = 'experiment/expmt_bottleneck24/out/source_pops_design.csv')
