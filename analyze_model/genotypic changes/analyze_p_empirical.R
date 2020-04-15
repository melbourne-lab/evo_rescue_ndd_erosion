# In this script: compare changes in allele frequencies from simulations with
# change in allele frequencies expected by some derived formulae.
# SN 10 Apr 2020

# Load in packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

# Clear your namespace
rm(list = ls())

# Source code needed to run simulations.
source('base_model_source/sim_functions.R')

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
                                   wfitn = sqrt(1 / 0.14/2),
                                   sig.e = 0))  # Selection strength

gt.l = test.sim %>%
  select(-c(r_i, fem)) %>%
  gather(key = allele.copy, value = val, -c(i, g_i, z_i, w_i, gen)) %>%
  mutate(locus = gsub('^[ab]', '', allele.copy)) %>%
  group_by(gen, locus, i) %>%
  mutate(n.pos = sum(val > 0)) %>%
  group_by(gen, locus, n.pos) %>%
  summarise(g.bar = mean(g_i),
            g.var = var(g_i),
            w.bar = mean(w_i),
            w.var = var(w_i),
            n = n()) %>%
  ungroup()

head(gt.l)

gt.l %>%
  mutate(n.jit = n.pos + rnorm(nrow(.), 0, 0.1)) %>%
  ggplot(aes(x = n.jit)) +
  geom_point(aes(y = g.bar, colour = gen),
             alpha = 0.75) +
  geom_segment(aes(xend = n.jit,
                   y = g.bar - 2 * sqrt(g.var / n),
                   yend = g.bar + 2 * sqrt(g.var / n)),
               alpha = 0.2) +
  geom_line(aes(y = g.bar, 
                group = interaction(gen, locus), 
                colour = gen),
            size = 0.125) +
  facet_wrap(~ locus)
  

gt.l %>%
  mutate(n.jit = n.pos + rnorm(nrow(.), 0, 0.1)) %>%
  ggplot(aes(x = n.jit)) +
  geom_segment(aes(x = -0.5, xend = 2.5,
                   y = 1, yend = 1),
               colour = 'gray55') +
  geom_point(aes(y = w.bar, colour = gen),
             alpha = 0.75) +
  geom_line(aes(y = w.bar, 
                group = interaction(gen, locus), 
                colour = gen),
            size = 0.125) +
  facet_wrap(~ locus)

gt.l %>%
  mutate(n.jit = n.pos + rnorm(nrow(.), 0, 0.1)) %>%
  ggplot(aes(x = n.jit)) +
  geom_segment(aes(x = -0.5, xend = 2.5,
                   y = 1, yend = 1),
               colour = 'gray55') +
  geom_point(aes(y = w.bar, fill = gen,
                 size = log(n)),
             shape = 21,
             alpha = 0.15) +
  geom_line(aes(y = w.bar, 
                group = interaction(gen, locus), 
                colour = gen),
            size = 0.125) +
  facet_wrap(~ locus)

gt.l %>%
  group_by(gen, locus) %>%
  mutate(w.norm = w.bar / sum(w.bar)) %>%
  ungroup() %>%
  mutate(n.jit = n.pos + rnorm(nrow(.), 0, 0.1)) %>%
  ggplot(aes(x = n.jit)) +
  geom_point(aes(y = w.norm, fill = gen,
                 size = log(n)),
             shape = 21,
             alpha = 0.15) +
  geom_line(aes(y = w.norm, 
                group = interaction(gen, locus), 
                colour = gen),
            size = 0.125) +
  facet_wrap(~ locus)
# Differences here are actually very small...

gt.t = test.sim %>%
  select(-c(r_i, fem)) %>%
  gather(key = allele.copy, value = val, -c(i, g_i, z_i, w_i, gen)) %>%
  mutate(locus = gsub('^[ab]', '', allele.copy)) %>%
  group_by(gen, locus, i) %>%
  mutate(n.pos = sum(val > 0)) %>%
  group_by(gen, locus) %>%
  mutate(p.pos = mean(val > 0)) %>%
  group_by(gen) %>%
  mutate(g.bar = mean(g_i),
         p.pos = mean(p.pos)) %>%
  group_by(gen, n.pos) %>%
  summarise(w.bar = mean(w_i),
            g.bar = g.bar[1],
            p = p.pos[1])

gt.t

gt.z = gt.t %>%
  mutate(m = 25) %>%
  mutate(n.pos = paste0('pos', n.pos)) %>%
  spread(key = n.pos, val = w.bar) %>%
  mutate(w_A = 2 * exp(-((g.bar - 2.6) + (1/(sqrt(25))))^2 / (2*sqrt(1/0.14/2)^2)),
         w_0 = 2 * exp(-((g.bar - 2.6))^2 / (2*sqrt(1/0.14/2)^2)),
         w_a = 2 * exp(-((g.bar - 2.6) - (1/(sqrt(25))))^2 / (2*sqrt(1/0.14/2)^2)))

gt.z %>%
  ggplot() +
  geom_point(aes(x = pos0, y = w_a, colour = gen), shape = 'triangle', size = 2) +
  geom_point(aes(x = pos1, y = w_0, colour = gen), shape = 'circle', size = 2) +
  geom_point(aes(x = pos2, y = w_A, colour = gen), shape = 'square', size = 2) +
  labs(x = 'Observed', y = 'Formula prediction')

gt.z %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = w_A), colour = 'red') +
  geom_line(aes(y = w_0), colour = 'blue') +
  geom_line(aes(y = w_a), colour = 'brown') +
  geom_line(aes(y = pos2), colour = 'red', linetype = 2) +
  geom_line(aes(y = pos1), colour = 'blue', linetype = 2) +
  geom_line(aes(y = pos0), colour = 'brown', linetype = 2) 

# but also, note wide range in predictions, tiny range in actual fitnesses...

gt.z %>%
  gather(key = n.pos, value = w.obs, -c(gen, g.bar, p, w_A, w_0, w_a)) %>%
  gather(key = coeff, value = w.prd, -c(gen, g.bar, p, n.pos, w.obs)) %>%
  filter((n.pos %in% 'pos0' & coeff %in% 'w_a') |
           (n.pos %in% 'pos1' & coeff %in% 'w_0') |
           (n.pos %in% 'pos2' & coeff %in% 'w_A'))

## Try looking at just mean number of offspring...

mt = test.sim %>%
  filter(fem) %>%
  select(-c(z_i, fem)) %>%
  gather(key = allele.copy, value = val, -c(i, g_i, w_i, r_i, gen)) %>%
  mutate(locus = gsub('^[ab]', '', allele.copy)) %>%
  group_by(gen, locus, i) %>%
  mutate(n.pos = sum(val > 0)) %>%
  group_by(gen, locus, n.pos) %>%
  summarise(g.bar = mean(g_i),
            g.var = var(g_i),
            w.bar = mean(w_i),
            w.var = var(w_i),
            r.bar = mean(r_i),
            n = n()) %>%
  ungroup()

head(mt)

mt %>%
  ggplot(aes(x = n.pos)) +
  geom_point(aes(y = r.bar, colour = gen, group = gen),
             position = position_dodge(width = 0.2)) +
  facet_wrap(~ locus)

mt %>%
  ggplot(aes(x = gen)) +
  geom_point(aes(y = r.bar, fill = factor(n.pos), group = n.pos),
             position = position_dodge(width = 0.3),
             size = 2, shape = 21) +
  scale_fill_manual(values = c('gray88', 'gray55', 'gray11')) +
  facet_wrap(~ locus) +
  theme(legend.position = 'none',
        panel.background = element_blank())

mt.cool = mt %>%
  select(gen, locus, n.pos, n) %>%
  mutate(n.pos = paste0('np', n.pos)) %>%
  spread(n.pos, n) %>%
  mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
  mutate(p.now = (np2 + 0.5 * np1) / (np0 + np1 + np2)) %>%
  merge(y = mt %>% 
          select(gen, locus, n.pos, r.bar) %>%
          group_by(gen, locus) %>%
          summarise(p.tp1 = sum(n.pos * r.bar / (2*sum(r.bar)))),
        by = c('gen', 'locus')) %>%
  merge(y = mt %>% 
          select(gen, locus, n.pos, r.bar) %>%
          mutate(n.pos = paste0('rp', n.pos)) %>%
          spread(n.pos, r.bar) %>%
          mutate_all(function(x) ifelse(is.na(x), 0, x)),
        by = c('gen', 'locus'))

mt.cool

mt.cool %>%
  ggplot() +
  geom_segment(x = 0, xend = 1, y = 0, yend = 1,
               colour = 'gray55') +
  geom_line(aes(x = p.now, y = p.tp1, group = locus)) +
  geom_point(aes(x = p.now, y = p.tp1, colour = gen)) +
  facet_wrap(~ locus) +
  theme(legend.position = 'none')


# I really have no idea what to do here.

# I really have no idea why this won't work as planned.
