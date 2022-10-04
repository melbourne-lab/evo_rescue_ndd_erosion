# Script for testing to see if genetic variance increases when allele
# frequencies are not initialized at parity.
# All else equal, genetic variance is maximized when p = 1/2.
# Here, we look to see if variance increases when starting with p < 1/2,
# with the optimum allele frequency at p > 1/2 such that variance should
# (in the absence of drift) increase during adaptation.
# This is meant for comparison with results of Barfield and Holt (2016)
# SN, 7 Apr 2022, rerun 25 Apr 2022, re-run and archived 3 Oct 2022

### Clear namespace
rm(list = ls())

### Load packages
library(ggplot2)
library(dplyr)
library(tidyr)

### Load source materials

# Get simulation functions
source('sim_functions.R')

# Define number of trials
trials = 1000

# Define parameters
pars = data.frame(
  trial = 1:trials,
  n.pop0 = 100,
  end.time = 16,
  init.row = 1e4,
  n.loci = 25,
  w.max = 2,
  theta = 1.4,
  wfitn = sqrt(3.5),
  sig.e = sqrt(0.5),
  pos.p = 0.36,
  alpha = 0
)

liszt = vector('list', nrow(pars))

set.seed(541002)

for (j in 1:nrow(pars)) {
  
  sim.output = sim( a = c(1/2, -1/2),
                    params = pars[j,])
  
  demo.summ = sim.output %>%
    group_by(gen) %>%
    summarise(n = n(),
              gbar = mean(g_i),
              zbar = mean(z_i),
              wbar = mean(w_i))
  
  gene.summ = sim.output %>%
    select(-c(g_i, z_i, w_i, r_i, fem)) %>%
    gather(key = loc.copy, value = val, -c(i, gen)) %>%
    mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
    group_by(i, gen, locus) %>%
    mutate(heter = !sum(val)) %>%
    group_by(gen, locus) %>%
    summarise(p = mean(val > 0),
              f = 1 - (mean(heter) / (2 * p * (1-p)))) %>%
    group_by(gen) %>%
    summarise(p.fix.pos = mean(p == 1),
              p.fix.neg = mean(p == 0),
              f = mean(f, na.rm = TRUE),
              v = sum(2 * p * (1 - p)) / pars$n.loci[1])
  
  liszt[[j]] = cbind(demo.summ, gene.summ %>% select(-gen)) %>% 
    mutate(trial = j)
  
  print(paste0(j, ' of ', nrow(pars)))
}

all.sims = do.call(rbind, liszt) %>%
  merge(y = pars %>% select(trial, n.pop0, alpha) %>% mutate(low.var = FALSE), by = 'trial') %>%
  group_by(trial)  %>%
  mutate(ext.gen = max(gen),
         extinct = !ext.gen %in% pars$end.time[1]) %>%
  ungroup()

all.sims %>%
  ggplot(aes(x = gen, y = v, group = trial)) +
  geom_line(aes(linetype = factor(extinct)), size = 0.1)

all.sims %>%
  ggplot(aes(x = gen, y = wbar, group = trial)) +
  geom_line(aes(linetype = factor(extinct)))

# Aggregated variance (v) and allele frequency (p)
vp = all.sims %>%
  # Re-index generation
  mutate(gen = gen - 1) %>%
  # Convert genotypes to allele frequencies
  mutate(p = 1/2 + gbar / (2 * sqrt(pars$n.loci[1]))) %>%
  group_by(gen) %>%
  summarise(
    vbar = mean(v),
    vvar = var(v),
    pbar = mean(p),
    pvar = var(p),
    nn   = n()
  ) #%>%

# 
vp %>%
  mutate(vopt = 2 * (pbar * (1 - pbar))) %>%
  ggplot(aes(x = gen)) +
  geom_ribbon(
    aes(
      ymin = vbar - 2*sqrt(vvar / nn),
      ymax = vbar + 2*sqrt(vvar / nn)
    ),
    alpha = 0.1
  ) +
  geom_segment(
    aes(
      x = 3.5, xend = 3.5,
      y = .4-.0125, yend = 0.5
    ),
    linetype = 2, colour = 'gray55'
  ) +
  geom_line(
    aes(
      y = vbar
    )
  ) +
  geom_line(
    aes(
      y = vopt,
    ), 
    colour = 'blue'
  ) +
  labs(x = 'Generation', y = 'Mean aditive genetic variance') +
  theme(panel.background = element_blank())

ggsave('analysis_results/fig_I1_incr_variance.png',
       width = 5, height = 5)

### Nothing below this line was included in final analysis

all.sims %>%
  ggplot(aes(x = gen, y = n, group = trial)) +
  geom_line() +
  scale_y_log10()

all.sims %>%
  ggplot(aes(x = gen, y = wbar, group = trial)) +
  geom_line() +
  scale_y_log10()

all.sims %>%
  ggplot(aes(x = gen, y = zbar, group = trial)) +
  geom_line()

all.sims %>%
  group_by(gen) %>%
  summarise(pos = mean(p.fix.pos),
            neg = mean(p.fix.neg),
            any.ext   = n() < 1000) %>%
  gather(key = allele, value = p, -c(gen, any.ext)) %>%
  ggplot(aes(x = gen, y = p, colour = allele)) +
  geom_line() +
  geom_point(aes(shape = any.ext), size = 3) +
  scale_shape_manual(values = c(19, 21)) +
  labs(y = 'Proportion of loci at fixation')
