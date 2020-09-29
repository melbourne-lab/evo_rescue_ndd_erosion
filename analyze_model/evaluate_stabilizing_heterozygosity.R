# Quick and dirty script to evalute changes in heterozygosity due to inbreeding
# load.
# Keller and Waller (2002) has a formula for deviation from expected
# heterozygosity, F = 1 - (h_obs)/(h_HWE).
# I want to see how easy this would be, and what the distribution of $F$, would
# look like under stabilizing selecton (i.e., no directional selection).
# SN September 28 2020

### Load packages
library(ggplot2)
library(dplyr)
library(tidyr)

### Clear namespace
rm(list = ls())

### Load source materials

# Get simulation functions
source('base_model_source/sim_functions.R')

# Define number of trials (per initial size)
trials = 20

# Define parameters
pars = expand.grid(local.trial = 1:trials,
                   n.pop0 = c(20, 60)) %>%
  mutate(end.time = 15,
         init.row = 1e4,
         n.loci = 25,
         w.max = 1.25,
         theta = 0,
         wfitn = sqrt(1 / 0.14 / 2),
         sig.e = sqrt(0.5),
         pos.p = 0.5)

### Run simulatons

liszt = vector('list', nrow(pars))

set.seed(4572899)

for (i in 1:nrow(pars)) {
  liszt[[i]] = sim( params = pars[i,])
  print(i)
}

# Put all sims into one object.
all.pops = unroller(liszt)

### First, take just one sim.

sim1 = liszt[[1]]

# Observed hertability
# (note: calculated by averaging loci over individuals, then averaging over
# loci)
sim1.obs.h = sim1 %>%
  select(-c(g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(i, gen)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(gen, locus, i) %>%
  summarise(h_i = !sum(val)) %>%
  group_by(gen, locus) %>%
  summarise(h_j = mean(h_i)) %>%
  group_by(gen) %>%
  summarise(h_t = mean(h_j))

# Heritabilty expected from HWE
# note: as with above, calculates heterzygosty per locus, then averages over
# loci.
sim1.hwe.h = sim1 %>%
  select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(gen)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(gen, locus) %>%
  summarise(p = mean(val > 0),
            h = 2 * p * (1 - p)) %>%
  group_by(gen) %>%
  summarise(h.pred = mean(h))

# Merge the two and calculate F: 1 - (h_observed) / (h_predicted)
sim1.inf = merge(sim1.obs.h, sim1.hwe.h,
                 by = c('gen')) %>%
  mutate(f = 1 - (h_t / h.pred))

# Plot f over time
sim1.inf %>%
  ggplot() +
  geom_point(aes(x = gen, y = f))

# Hmm... not what I was expecting. All negative.
# Of course this is just looking at deviation from random mating.
# Semems like there's slightly _more_ heterozygosity than expected.

### Try all sims.

# Observed heterozygosity (same as above).
all.obs.h = all.pops %>%
  select(-c(g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(i, gen, trial)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(trial, gen, locus, i) %>%
  summarise(h_i = !sum(val)) %>%
  group_by(trial, gen, locus) %>%
  summarise(h_j = mean(h_i)) %>%
  group_by(trial, gen) %>%
  summarise(h_t = mean(h_j))

# Expected heterozygosity (same as above)
all.hwe.h = all.pops %>%
  select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(gen, trial)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(trial, gen, locus) %>%
  summarise(p = mean(val > 0),
            h = 2 * p * (1 - p)) %>%
  group_by(trial, gen) %>%
  summarise(h.pred = mean(h))

# Merge
all.inf = merge(all.obs.h, all.hwe.h,
                 by = c('trial', 'gen')) %>%
  mutate(f = 1 - (h_t / h.pred))

# Plot of all trials
# Black = survving to end, red = extinct
all.inf %>%
  group_by(trial) %>% mutate(surv = max(gen) == 15) %>% ungroup() %>%
  ggplot() +
  geom_line(aes(x = gen, y = f, group = trial, colour = surv)) +
  scale_color_manual(values = c('red', 'black'))

# Quick sketch: are these different from zero?
all.inf %>%
  group_by(gen) %>%
  summarise(mu = mean(f),
            sig = var(f),
            mup2 = mu + 2 * sqrt(sig / 40))
# All negative, all look different from zero at 95% level

# Conclusion: definitely no evidence to conclude that heterozygosity is _less_
# than expected. This result is not totally surprising as mating here is totally
# random. It's possible there's slightly _more_ heterozygosity than expected.
