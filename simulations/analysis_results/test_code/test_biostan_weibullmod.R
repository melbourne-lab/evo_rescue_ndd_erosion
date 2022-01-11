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

### Looking at other code snippets from tutorial
# (https://rstudio-pubs-static.s3.amazonaws.com/435225_07b4ab5afa824342a4680c9fb2de6098.html)

bayesplot::mcmc_acf(as.data.frame(stan2.fit), 
                    pars = c('mu', 'alpha', 'beta_bg[1]', 'beta_bg[2]', 'beta_bg[3]',
                             'beta_bg[4]', 'beta_bg[5]', 'beta_bg[6]', 'beta_bg[7]'))
# This actually looks good.

# Oh the remaining stuff is hard coded into the stan file
# Really don't want to run this with n in the thousands and produce thousands of columns

# I suppose we can generate survival curves just from these posterior draws
# (n.b. they don't include prediction error, i.e., are just estimates of the mean)

stan2.wide = as.data.frame(stan2.fit) %>%
  rename(beta1 = `beta_bg[1]`,
         beta2 = `beta_bg[2]`,
         beta3 = `beta_bg[3]`,
         beta4 = `beta_bg[4]`,
         beta5 = `beta_bg[5]`,
         beta6 = `beta_bg[6]`,
         beta7 = `beta_bg[7]`) %>%
  select(alpha, mu, contains('beta')) %>%
  select(-contains('bg'))

post.draws = merge(
  x = stan2.wide %>% mutate(i = 1:nrow(.)) %>% filter(i < 201),
  y = data.frame(expand.grid(i = 1:200, t = 1:50))
) %>%
  arrange(i, t) %>%
  mutate(sig000 = exp(-mu / alpha),
         sig100 = exp(-(mu + beta1) / alpha),
         sig010 = exp(-(mu + beta2) / alpha),
         sig110 = exp(-(mu + beta1 + beta2 + beta4) / alpha),
         sig001 = exp(-(mu + beta3) / alpha),
         sig101 = exp(-(mu + beta1 + beta3 + beta6) / alpha),
         sig011 = exp(-(mu + beta2 + beta3 + beta5) / alpha),
         sig111 = exp(-(mu + beta1 + beta2 + beta3 + 
                        beta4 + beta5 + beta6 + beta7) / alpha),
         st_a0v0s0 = exp(-t / sig000)^alpha,
         st_a1v0s0 = exp(-t / sig100)^alpha,
         st_a0v1s0 = exp(-t / sig010)^alpha,
         st_a1v1s0 = exp(-t / sig110)^alpha,
         st_a0v0s1 = exp(-t / sig001)^alpha,
         st_a1v0s1 = exp(-t / sig101)^alpha,
         st_a0v1s1 = exp(-t / sig011)^alpha,
         st_a1v1s1 = exp(-t / sig111)^alpha) %>%
  select(i, t, contains('st_')) %>%
  gather(key = st, value = estim, -c(i, t)) %>%
  mutate(ndd = grepl('a1', st),
         lov = grepl('v1', st),
         lrg = grepl('s1', st))

post.draws %>%
  ggplot(aes(x = t, y = estim)) +
  geom_line(
    aes(
        group = interaction(i, st), 
        colour = ndd
      ),
    size = 0.1
    ) +
  scale_colour_manual(values = c('black', 'purple')) +
  facet_wrap(lov ~ lrg)

# Cool, intuitive, etc.

# Although, what we are actually interested in is the mean time until extinction.

post.preds = stan2.wide %>% 
  mutate(i = 1:nrow(.)) %>% 
  mutate(siga0v0s0 = exp(-mu / alpha),
         siga1v0s0 = exp(-(mu + beta1) / alpha),
         siga0v1s0 = exp(-(mu + beta2) / alpha),
         siga1v1s0 = exp(-(mu + beta1 + beta2 + beta4) / alpha),
         siga0v0s1 = exp(-(mu + beta3) / alpha),
         siga1v0s1 = exp(-(mu + beta1 + beta3 + beta6) / alpha),
         siga0v1s1 = exp(-(mu + beta2 + beta3 + beta5) / alpha),
         siga1v1s1 = exp(-(mu + beta1 + beta2 + beta3 + 
                          beta4 + beta5 + beta6 + beta7) / alpha)) %>%
  select(-contains('beta')) %>%
  select(-mu) %>%
  gather(key = sig, value = sig.estim, -c(i, alpha)) %>%
  mutate(ndd = grepl('a1', sig),
         lov = grepl('v1', sig),
         lrg = grepl('s1', sig)) %>%
  ungroup() %>%
  mutate(pred = rweibull(n = nrow(.), shape = alpha, scale = sig.estim))

post.preds %>%
  ggplot(aes(x = pred)) +
  geom_density(
    aes(
      colour = ndd
    )
  ) +
  scale_x_log10(limits = c(1, 50)) +
  facet_wrap(lov ~ lrg)
# Not a fan of this! Wrong denominator when censoring.  

# Are mean times to extinction relevant here? Growth rates increasing...
# Can we extrapolate beyond 50?

# Let's actually compare our model with actual data.

# Empirical survival curves probably look like this?
# (note: can't do this with raw data due to size censoring)
#p.surv = 
p.surv = ext.times %>% 
  filter(extinct) %>% 
  group_by(n.pop0, low.var, alpha, ext.gen) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  rbind(
    expand.grid(
      n.pop0 = c(20, 100),
      low.var = c(TRUE, FALSE),
      alpha = c(0, 0.0035),
      ext.gen = 1:50
    ) %>%
    mutate(n = 0)
  ) %>%
  group_by(n.pop0, low.var, alpha, ext.gen) %>%
  summarise(n = sum(n)) %>%
  group_by(n.pop0, low.var, alpha) %>%
  mutate(n.ext = cumsum(n)) %>% 
  ungroup() %>%
  mutate(p.surv = 1 - n.ext/1000)
  
p.surv %>%
  ggplot(aes(x = ext.gen, y = p.surv)) +
  geom_line(
    aes(
      colour = factor(alpha)
    ),
  ) +
  scale_colour_manual(values = c('black', 'purple')) +
  facet_wrap(low.var ~ n.pop0)

# Oh... duh I have plots like this already
# Anyway this doesn't really look like the model estimates
# the hazard rate really needs to fall a lot more quickly
# (different alpha value)
