# What happens when analyzing extinction times with permutation testing?
# Not sure if censoring (and therefore different sample sizes) will bias things?

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

all.data = 
  rbind(
    read.csv('simulations/outputs/longsims/longsims_n100_a000_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n100_a000_lowvar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a000_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a000_lowvar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n100_a035_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n100_a035_lowvar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a035_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a035_lowvar.csv')
  )

nrow(all.data)

# Add extinction column
all.data = all.data %>%
  arrange(trial, n.pop0, low.var, alpha, desc(gen)) %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  group_by(n.pop0, low.var, alpha, trial) %>%
  summarise(extinct = gen < 50 & n < 10000) %>%
  merge(x = all.data, y = ., by = c('trial', 'n.pop0', 'low.var', 'alpha'))


# Get just data on extinctions.
ext.times = all.data %>%
  # Get time of extinction
  group_by(trial, n.pop0, low.var, alpha) %>%
  summarise(ext.gen = max(gen),
            extinct = extinct[1]) %>%
  # Do this because of censoring of population size
  # (I stopped simulations when populations grew too large under assumption they would not go extinct)
  mutate(ext.gen = ifelse(ext.gen < 50 & !extinct, 50, ext.gen)) %>%
  filter(extinct)

# Observed differences in extinction time between NDD and DI populations
obs.diffs = ext.times %>%
  group_by(n.pop0, low.var) %>%
  summarise(ext.diff = mean(ext.gen[alpha > 0]) - mean(ext.gen[alpha == 0])) %>%
  ungroup()

nul.diffs = vector(mode = 'list')# expand.grid(i = 1:1000, n.pop0 = c(20, 100), low.var = c(TRUE, FALSE))

for (trial in 1:1000) {
  nul.diffs[[trial]] = ext.times %>%
    group_by(n.pop0, low.var) %>%
    mutate(ext.gen = sample(ext.gen)) %>%
    summarise(ndiff = mean(ext.gen[alpha > 0]) - mean(ext.gen[alpha == 0])) %>%
    ungroup()
}

merge(
  x = do.call(nul.diffs, what = rbind) %>% group_by(n.pop0, low.var),
  y = obs.diffs
) %>%
  group_by(n.pop0, low.var) %>%
  summarise(p = mean(ext.diff > ndiff))

# wait... we're probably more interested in CIs than a p-value

ci.diffs = vector(mode = 'list')# expand.grid(i = 1:1000, n.pop0 = c(20, 100), low.var = c(TRUE, FALSE))

for (trial in 1:1000) {
  ci.diffs[[trial]] = ext.times %>%
    group_by(n.pop0, low.var, alpha) %>%
    mutate(ext.gen = sample(ext.gen, replace = TRUE)) %>%
    group_by(n.pop0, low.var) %>%
    summarise(sdiff = mean(ext.gen[alpha > 0]) - mean(ext.gen[alpha == 0])) %>%
    ungroup()
}

do.call(ci.diffs, what = rbind) %>% 
  group_by(n.pop0, low.var) %>%
  summarise(ci025 = quantile(sdiff, 0.025),
            ci975 = quantile(sdiff, 0.975))

