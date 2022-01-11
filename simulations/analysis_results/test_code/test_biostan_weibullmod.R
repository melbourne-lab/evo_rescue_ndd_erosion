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
                 seed = 10,
                 cores = 4)

stan0.fit 
# this is better
# does very poorly for large popupations...  my guess is because there's a lot of censoring

# taus: I don't think these are interesting...
# what are we actually interested in? I guess alpha, mu, beta_bg?

stan0.pars = as.data.frame(stan0.fit) %>%
  rename(beta = `beta_bg[1]`) %>%
  mutate(i = rep(1:4000)) %>%
  gather(key = par, value = estim, -i) %>%
  mutate(chain = (i-1) %/% 1000 + 1,
         iter  = (i-1) %% 1000 + 1)
  
ggplot(stan0.pars %>% filter(par %in% c('mu', 'alpha', 'beta'))) +
  geom_line(aes(x = iter, y = estim, group = chain),
            size = 0.5) +
  facet_wrap(chain ~ par, scales = 'free_y', nrow = 4)

# Chains 1 and 3 look bad for a few iterations

bayesplot::mcmc_areas(stan0.fit, pars = c('mu', 'alpha', 'beta_bg[1]'))

### Now, can this code work with a two-way interaction?

ext1 = ext.times %>% filter(n.pop0 < 50)

mod.mat1 = ext1 %>%
  ungroup() %>%
  mutate(b_alpha = as.numeric(alpha > 0),
         b_lovar = as.numeric(low.var),
         b_alpha_lovar = b_alpha * b_lovar) %>%
  select(b_alpha, b_lovar, b_alpha_lovar)

stan.inputs1 = list(
  Nobs = sum(ext1$extinct),
  Ncen = sum(!ext1$extinct),
  M_bg = 3,
  yobs = ext1$ext.gen[ext1$extinct],
  ycen = ext1$ext.gen[!ext1$extinct],
  Xobs_bg = as.matrix(mod.mat1[ext1$extinct,]),
  Xcen_bg = as.matrix(mod.mat1[!ext1$extinct,])
)

stan1.fit = stan(file = 'simulations/weibull_survival_model.stan',
                 data = stan.inputs1,
                 seed = 10,
                 cores = 4)

# The code works

stan1.fit
# what to conclude when ESS for a raw parameter is low but is high for transformed??

stan1.pars = as.data.frame(stan1.fit) %>%
  rename(beta1 = `beta_bg[1]`,
         beta2 = `beta_bg[2]`,
         beta3 = `beta_bg[3]`) %>%
  mutate(i = rep(1:4000)) %>%
  gather(key = par, value = estim, -i) %>%
  mutate(chain = (i-1) %/% 1000 + 1,
         iter  = (i-1) %% 1000 + 1)

ggplot(stan1.pars %>% filter(par %in% c('mu', 'alpha', 'beta1', 'beta2', 'beta3'))) +
  geom_line(aes(x = iter, y = estim, group = chain),
            size = 0.5) +
  facet_wrap(chain ~ par, scales = 'free_y', nrow = 4)
# wtf is happening in chain 4...?

bayesplot::mcmc_areas(as.data.frame(stan1.fit)[1:2800,], 
                      pars = c('mu', 'alpha', 'beta_bg[1]', 'beta_bg[2]', 'beta_bg[3]'))

### Try a three-way model??

mod.mat2 = ext.times %>%
  ungroup() %>%
  mutate(b_alpha = as.numeric(alpha > 0),
         b_lovar = as.numeric(low.var),
         b_large = as.numeric(n.pop0 > 50),
         b_alpha_lovar = b_alpha * b_lovar,
         b_lovar_large = b_lovar * b_large,
         b_alpha_large = b_alpha * b_large,
         b_alpha_lovar_large = b_alpha * b_lovar * b_large) %>%
  select(contains('b_'))

stan.inputs2 = list(
  Nobs = sum(ext.times$extinct),
  Ncen = sum(!ext.times$extinct),
  M_bg = 7,
  yobs = ext.times$ext.gen[ext.times$extinct],
  ycen = ext.times$ext.gen[!ext.times$extinct],
  Xobs_bg = as.matrix(mod.mat2[ext.times$extinct,]),
  Xcen_bg = as.matrix(mod.mat2[!ext.times$extinct,])
)

stan2.fit = stan(file = 'simulations/weibull_survival_model.stan',
                 data = stan.inputs2,
                 seed = 10,
                 cores = 4)

stan2.fit
# whoa... all of these converged??
# is censoring dependent on parameters?

stan2.pars = as.data.frame(stan2.fit) %>%
  rename(beta1 = `beta_bg[1]`,
         beta2 = `beta_bg[2]`,
         beta3 = `beta_bg[3]`,
         beta4 = `beta_bg[4]`,
         beta5 = `beta_bg[5]`,
         beta6 = `beta_bg[6]`,
         beta7 = `beta_bg[7]`) %>%
  mutate(i = rep(1:4000)) %>%
  gather(key = par, value = estim, -i) %>%
  mutate(chain = (i-1) %/% 1000 + 1,
         iter  = (i-1) %% 1000 + 1)

ggplot(stan2.pars %>% filter(par %in% c('mu', 'alpha', 'beta1', 'beta2', 'beta3',
                                        'beta4', 'beta5', 'beta6', 'beta7'))) +
  geom_line(aes(x = iter, y = estim, group = chain),
            size = 0.5) +
  facet_wrap(par ~ chain, scales = 'free_y', ncol = 4)

bayesplot::mcmc_areas(as.data.frame(stan2.fit), 
                      pars = c('mu', 'alpha', 'beta_bg[1]', 'beta_bg[2]', 'beta_bg[3]',
                               'beta_bg[4]', 'beta_bg[5]', 'beta_bg[6]', 'beta_bg[7]'))
# certainties... maybe this is fine

