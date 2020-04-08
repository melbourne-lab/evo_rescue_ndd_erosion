# In this script: compare changes in allele frequencies from simulations with
# change in allele frequencies expected by some derived formulae.
# SN 7 Apr 2020

# Load in packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

# Clear your namespace
rm(list = ls())

# Source code needed to run simulations.
source('init_model/sim_evolve_allele_spec_brs.R')

# Set seed for reproducibility
set.seed(40904)

# Run one simulation
test.sim = sim(a = c(-1/2, 1/2),
               params = data.frame(end.time = 15,      # Length of simulatin
                                   init.row = 1e4,     # Number of rows to initialize data
                                   n.loci = 25,        # Number of loci
                                   n.pop0 = 20,        # Initial population size
                                   w.max = 2,          # Max fitness of an individual
                                   theta = 2.6,        # Distance from phenotypic optimum
                                   w.sqd = 1 / 0.14))  # Selection strength

gt.t = test.sim %>%
  select(-c(r_i, fem)) %>%
  gather(key = allele.copy, value = val, -c(i, b_i, w_i, gen)) %>%
  mutate(locus = gsub('^[ab]', '', allele.copy)) %>%
  group_by(gen, locus, i) %>%
  mutate(n.pos = sum(val > 0)) %>%
  group_by(gen, locus) %>%
  mutate(p.pos = mean(val > 0)) %>%
  group_by(gen) %>%
  mutate(b.bar = mean(b_i),
         p.pos = mean(p.pos)) %>%
  group_by(gen, n.pos) %>%
  summarise(w.bar = mean(w_i),
            b.bar = b.bar[1],
            p = p.pos[1])

gt.t

gt.z = gt.t %>%
  mutate(n.pos = paste0('pos', n.pos)) %>%
  spread(key = n.pos, val = w.bar) %>%
  mutate(w_A = 2 * exp(-((b.bar - 2.6) + (1/sqrt(25)))^2 / (1/0.14)),
         w_0 = 2 * exp(-((b.bar - 2.6))^2 / (1/0.14)),
         w_a = 2 * exp(-((b.bar - 2.6) - (1/sqrt(25)))^2 / (1/0.14)))

gt.z %>%
  ggplot() +
  geom_point(aes(x = pos0, y = w_a, colour = gen), shape = 'triangle', size = 2) +
  geom_point(aes(x = pos1, y = w_0, colour = gen), shape = 'circle', size = 2) +
  geom_point(aes(x = pos2, y = w_A, colour = gen), shape = 'square', size = 2) +
  labs(x = 'Observed', y = 'Formula prediction')

gt.z %>%
  ggplot() +
  geom_point(aes(x = pos0, y = w_a, colour = gen), shape = 'triangle', size = 2) +
  geom_point(aes(x = pos1, y = w_0, colour = gen), shape = 'circle', size = 2) +
  geom_point(aes(x = pos2, y = w_A, colour = gen), shape = 'square', size = 2) +
  geom_line(aes(x = pos0, y = w_a)) +
  geom_line(aes(x = pos1, y = w_0)) +
  geom_line(aes(x = pos2, y = w_A)) +
  labs(x = 'Observed', y = 'Formula prediction')

gt.z %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = w_A), colour = 'red') +
  geom_line(aes(y = w_0), colour = 'blue') +
  geom_line(aes(y = w_a), colour = 'brown') +
  geom_line(aes(y = pos2), colour = 'red', linetype = 2) +
  geom_line(aes(y = pos1), colour = 'blue', linetype = 2) +
  geom_line(aes(y = pos0), colour = 'brown', linetype = 2) 

# wtf is happening in last time step...?
# but also, note wide range in predictions, tiny range in actual fitnesses...

gt.z %>%
  gather(key = n.pos, value = w.obs, -c(gen, b.bar, p, w_A, w_0, w_a)) %>%
  gather(key = coeff, value = w.prd, -c(gen, b.bar, p, n.pos, w.obs)) %>%
  filter((n.pos %in% 'pos0' & coeff %in% 'w_a') |
         (n.pos %in% 'pos1' & coeff %in% 'w_0') |
         (n.pos %in% 'pos2' & coeff %in% 'w_A'))
  ggplot() +
  geom_line(aes(x = ))