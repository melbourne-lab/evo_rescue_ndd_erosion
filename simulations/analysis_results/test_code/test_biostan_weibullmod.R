# Sandbox script for looking at biostan weibull model
# SN 10 Jan 2022
# biostan script: https://github.com/jburos/biostan/blob/master/inst/stan/weibull_survival_model.stan

library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)

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

# Actually... this may be sufficient.for testing the weibull code

# Get just data on extinctions.
ext.times = all.data %>%
  # Get time of extinction
  group_by(trial, n.pop0, low.var, alpha) %>%
  summarise(ext.gen = max(gen),
            extinct = extinct[1]) %>%
  # Do this because of censoring of population size
  # (I stopped simulations when populations grew too large under assumption they would not go extinct)
  mutate(ext.gen = ifelse(ext.gen < 50 & !extinct, 50, ext.gen))

with(ext.times, table(extinct, ext.gen))

# Now, we only want to look at one parameter.
# Say, alpha.
# Get only one parameter combination for remaining params

ext0 = ext.times %>% filter(n.pop0 < 50, low.var)

nrow(ext0)
ext0 %>% head()

# Okay. Censoring is in the extinct column

### Try the stan model

stan.inputs0 = list(
  Nobs = sum(ext0$extinct),
  Ncen = sum(!ext0$extinct),
  M_bg = 1,
  yobs = ext0$ext.gen[ext0$extinct],
  ycen = ext0$ext.gen[!ext0$extinct],
  Xobs_bg = matrix(as.numeric(ext0$alpha[ext0$extinct] > 0)),
  Xcen_bg = matrix(as.numeric(ext0$alpha[!ext0$extinct] > 0))
)

stan0.fit = stan(file = 'simulations/weibull_survival_model.stan',
                 data = stan.inputs0,
                 cores = 4)

stan0.fit 
# this is better
# does very poorly for large popupations...  my guess is because there's a lot of censoring

# taus: I don't think these are interesting...
# what are we actually interested in? I guess alpha, mu, beta_bg?