# A script to subsample a uniform number of individuals out of each of the
# bottlenecked populations.
# Doing this to ensure that there is no among-source population differences in
# sampling variation (due to differences in source population size).
# SN - August 28 2020

# I'm using 600 as a subsample size.
# It's neatly divisible by both 60 and 20.

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Load in all 48 thrice-bottlenecked source populations
neck.popns = read.csv('simulations/prepare_simulations/outputs/thrice_bottlenecked_popns.csv',
                      stringsAsFactors = FALSE)

# First: remove all populations that are smaller than 600.
# (Can't subsample 600 individuals from a population smaller than 600)
cens.necks = neck.popns %>%
  group_by(trial) %>%
  filter(n() > 600) %>%
  ungroup()

# Did this work?
head(cens.necks)

# There are now 27 source populations.
unique(cens.necks$trial) %>% length()

# Take a look at the genetic variation in these populations.
cens.vars = cens.necks %>%
  select(-c(i, g_i, z_i, w_i, r_i, fem, gen)) %>%
  gather(key = loc.copy, value = val, -c(trial)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(trial, locus) %>%
  summarise(p = mean(val > 0)) %>%
  group_by(trial) %>%
  summarise(p.fix.pos = mean(p == 1),
            p.fix.neg = mean(p == 0),
            v = sum(2 * p * (1 - p)) / pars$n.loci[1]) %>%
  group_by(trial) %>%
  summarise(p.fix.pos = mean(p.fix.pos),
            p.fix.neg = mean(p.fix.neg),
            vbar = mean(v))
# Plot
ggplot(cens.vars) + geom_histogram(aes(x = vbar), binwidth = 0.025)

# Surprisingly, still kind of low.
mean(cens.vars$vbar)

# Now, do the subsampling. Take 600 individuals from each population.

# Set seed
set.seed(25025)

# Do subsampling
cens.samps = cens.necks %>%
  group_by(trial) %>%
  sample_n(size = 600) %>%
  ungroup()

# Write out data.
write.csv(cens.samps, row.names = FALSE,
          file = 'simulations/prepare_simulations/outputs/thrice_bottlenecked_600s.csv')

# Inspect the variance in these populations (did subampling affect variance?)
cens.samps %>%
  select(-c(i, g_i, z_i, w_i, r_i, fem, gen)) %>%
  gather(key = loc.copy, value = val, -c(trial)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(trial, locus) %>%
  summarise(p = mean(val > 0)) %>%
  group_by(trial) %>%
  summarise(p.fix.pos = mean(p == 1),
            p.fix.neg = mean(p == 0),
            v = sum(2 * p * (1 - p)) / pars$n.loci[1]) %>%
  group_by(trial) %>%
  summarise(p.fix.pos = mean(p.fix.pos),
            p.fix.neg = mean(p.fix.neg),
            vbar = mean(v)) %>%
  ggplot() +
  geom_histogram(aes(x = vbar), binwidth = 0.025)

# Not visibly.
