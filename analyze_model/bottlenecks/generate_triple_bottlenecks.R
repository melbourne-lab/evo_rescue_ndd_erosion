##### Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

##### Clear namespace
rm(list = ls())

##### Load in source functions
source('base_model_source/sim_functions.R')

##### Parameters

# I'll run 2000 trials here.

# Use 25 loci for now (initial trials suggest this is large enough)
# theta = 0 - populations here are already adapted to their population
# wfitn is the same value used in other situations


n.trials = 200
n.necks = 3
thresh = 100

n.loci = 25
n.gens = 10
theta = 0
wfitn = sqrt(1 / 0.14 / 2)

params = data.frame(end.time = n.gens,
                    init.row = 1e5,
                    n.pop0 = 2,
                    n.loci = n.loci,
                    w.max = 2,
                    theta = theta,
                    wfitn = wfitn,
                    alpha = 0,
                    sig.e = 0)

liszt = vector('list', n.trials)

for (trial in 1:n.trials) {
  #print(paste0('trial', trial))
  base.pop = init.sim(params = params) %>%
    mutate(fem = c(TRUE, FALSE),
           r_i = ifelse(fem, rpois(n = 1, lambda = 2 * w_i), 0))
  if (sum(base.pop$r_i) > 1) {
    for (neck in 1:n.necks) {
      if (nrow(base.pop) > 1) {
        prop.pop = propagate.sim(params = params,
                                 popn = base.pop,
                                 evolve = TRUE)
        if (nrow(prop.pop) > 1) {
          base.pop =  prop.pop %>%
            sample_n(size = 2, replace = FALSE) %>%
            mutate(fem = c(TRUE, FALSE))
          #print(paste0('bottleneck', neck))
        } else { break }
      }
    }
  }
  if (neck == n.necks) {
    liszt[[trial]] = sim(
      params = params,
      init.popn = base.pop
    )
  }
  print(trial)
}

unrolled = unroller(liszt)

table(unrolled$trial)

demo.stats = unrolled %>%
  group_by(trial, gen) %>%
  summarise(n = n(),
            gbar = mean(g_i),
            gvar = var(g_i),
            wbar = mean(w_i),
            wvar = var(w_i)) %>%
  group_by(trial) %>%
  mutate(extinct = max(gen) < n.gens) %>%
  ungroup()

gene.stats = unrolled %>%
  select(-c(g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = locus.copy, value = gval, -c(i, gen, trial)) %>%
  mutate(locus = gsub('^[ab]', '', locus.copy)) %>%
  group_by(trial, gen, locus) %>%
  summarise(p = mean(gval > 0)) %>%
  group_by(trial, gen) %>%
  summarise(v = 2 * sum(p * (1-p)) / n.loci) %>%
  group_by(trial) %>%
  mutate(extinct = max(gen) < n.gens) %>%
  ungroup()
  

demo.stats %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = n,
                group = trial,
                colour = extinct)) +
  scale_y_log10() +
  scale_color_manual(values = c('black', 'red'))

gene.stats %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = v,
                group = trial,
                colour = extinct)) +
  scale_color_manual(values = c('black', 'red'))

gene.stats %>%
  group_by(gen, extinct) %>%
  mutate(vbar.ext = mean(v),
         vvar.ext = var(v),
         n.by.ext = n()) %>%
  group_by(gen) %>%
  mutate(vbar = mean(v),
         vvar = var(v),
         n.trials = n()) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = v,
                group = trial,
                colour = extinct)) +
  geom_line(aes(y = vbar.ext,
                group = extinct,
                colour = extinct),
            size = 4,
            linetype = 2) +
  geom_ribbon(aes(ymin = vbar.ext - 2 * sqrt(vvar.ext/n.by.ext),
                  ymax = vbar.ext + 2 * sqrt(vvar.ext/n.by.ext),
                  group = extinct,
                  fill = extinct),
              alpha = 0.2) +
  scale_color_manual(values = c('black', 'red')) +
  scale_fill_manual(values = c('black', 'red'))

demo.stats %>% filter(gen %in% n.gens) %>%
  ggplot() + geom_histogram(aes(x = n))

unrolled %>%
  filter(gen %in% n.gens) %>%
  write.csv(file = 'analyze_model/bottlenecks/thrice_bottlenecked_popns.csv',
            row.names = FALSE)
