# More code seeing why the rate of phenotypic change doesn't match Lande.

### Clear namespace
rm(list = ls())

### Load packages

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

### Load source materials

# Get simulation functions
source('base_model_source/sim_functions.R')

# Define number of trials
trials = 4000

# Define parameters
pars = expand.grid(local.trial = 1:trials) %>%
  mutate(global.trial = 1:nrow(.),
         end.time = 2,
         init.row = 300,
         n.loci = 25,
         n.pop0 = 60,
         w.max = 2,
         theta = 2.75,
         wfitn = sqrt(1 / 0.14 / 2),
         sig.e = sqrt(0.5),
         pos.p = 0.5,
         alpha = 0)

liszt = vector('list', nrow(pars))

set.seed(72)

for (i in 1:nrow(pars)) {
  sim.output = sim( a = c(1/2, -1/2), params = pars[i,])
  
  demo.summ = sim.output %>%
    group_by(gen) %>%
    summarise(n = n(),
              gbar = mean(g_i),
              gvar = var(g_i),
              zbar = mean(z_i),
              zvar = var(z_i),
              evar = var(z_i - g_i))
  
  gene.summ = sim.output %>%
    select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
    gather(key = loc.copy, value = val, -gen) %>%
    mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
    group_by(gen, locus) %>%
    summarise(p = mean(val > 0)) %>%
    group_by(gen) %>%
    summarise(v = sum(2 * p * (1 - p)) / pars$n.loci[1])
  
  liszt[[i]] = cbind(demo.summ, gene.summ %>% select(-gen)) 
  
  print(paste0(i, ' of ', nrow(pars)))
}

all.pops = unroller(liszt)

all.pops %>%
  ggplot() +
  geom_point(aes(x = v, y = gvar))

all.pops %>%
  filter(gen %in% 1) %>%
  ggplot() +
  geom_point(aes(x = zvar, y = evar + gvar)) +
  geom_segment(aes(x = 0.5, xend = 1.5, y = 0.5, yend = 1.5),
               colour = 'blue')

# aha... something hither is fucked

lm(zvar ~ sumvar, all.pops %>% filter(gen < 2) %>% mutate(sumvar = evar + gvar)) %>% summary()

# wait shit I've done this before
# ahhh what teh fuck

lm(zvar ~ sumvar, all.pops %>% filter(gen < 2) %>% mutate(sumvar = evar + gvar)) %>% plot()

all.pops %>%
  filter(gen %in% 1) %>%
  ggplot(aes(x = evar + gvar, y = zvar)) +
  geom_point() +
  geom_segment(aes(x = 0.5, xend = 1.5, y = 0.5, yend = 1.5),
               colour = 'orange') +
  stat_smooth(method = 'lm')

# oh mannnn was I just.. doing the wrong axes? no seriously how the fuck did that happen?
# I have never beeen so confused

all.pops %>%
  ggplot() +
  geom_point(aes(x = gvar, y = evar)) +
  geom_segment(aes(x = 0.3, xend = 1.0, y = 0.3, yend = 1.0),
               colour = 'blue')
# this is a little funky but looks fine

all.pops %>%
  filter(gen < 2) %>%
  ggplot() +
  geom_histogram(aes(x = evar), binwidth = 0.01) +
  geom_segment(aes(x = 0.5, xend = 0.5, y = 0, yend = 200),
               colour = 'blue')

# Yep, that looks good.

### Lande comparison

w2 = pars$wfitn

comp = all.pops %>%
  mutate(th = pars$theta[1]) %>%
  group_by(trial) %>%
  mutate(d1 = th - zbar[1],
         d2 = th - zbar[2]) %>%
  ungroup() %>%
  filter(gen < 2)

runner = comp %>%
  mutate(lande1 = d1 * (pars$wfitn[1]^2 + evar) / (pars$wfitn[1]^2 + evar + gvar),
         lande2 = d1 * (pars$wfitn[1]^2 + evar) / (pars$wfitn[1]^2 + evar + v),
         h2 = gvar / zvar,
         h2_vz = v / zvar,
         h2_ve = v / (v + evar),
         h2_ge = gvar / (gvar + evar),
         lande3 = d1 * (1 - h2 * evar / (evar + pars$wfitn[1]^2)), # note: mistake in here...
         lande4 = d1 * (1 - h2 * zvar / (zvar + pars$wfitn[1]^2)),
         lande5 = d1 * (1 - h2_ge * zvar / (zvar + pars$wfitn[1]^2)),
         lande6 = d1 * (1 - h2_vz * zvar / (zvar + pars$wfitn[1]^2)),
         lande7 = d1 * (pars$wfitn[1]^2 + evar) / (pars$wfitn[1]^2 + zvar),
         lande8 = d1 * (1 - (1-h2) * zvar / (zvar + pars$wfitn[1]^2)))

ggplot(comp, aes(x = lande1, y = d2)) + 
  geom_point() +
  geom_segment(aes(x = 2, xend = 3, y = 2, yend = 3), colour = 'blue')

# okay. What to do here.
# the lande prediction is consistently below 

ggplot(runner, aes(x = lande2, y = d2)) + 
  geom_point() +
  geom_segment(aes(x = 2, xend = 3, y = 2, yend = 3), colour = 'blue')

ggplot(comp, aes(x = h2_v, y = h2)) + 
  geom_point() +
  geom_segment(aes(x = 0.25, xend = 1, y = 0.25, yend = 1), colour = 'blue')
# no evidence here that h is off

ggplot(comp, aes(x = lande3, y = d2)) + 
  geom_point() +
  geom_segment(aes(x = 2, xend = 3, y = 2, yend = 3), colour = 'blue')
# aha. This one is good.
# wait... but it shouldn't be? there's a mistake in this one...

ggplot(comp, aes(x = lande4, y = d2)) + 
  geom_point() +
  geom_segment(aes(x = 2, xend = 3, y = 2, yend = 3), colour = 'blue')
# uhm?

# wtf?

ggplot(runner, aes(x = lande3 - d2)) + geom_histogram(, binwidth = 0.01)
# (still somewhat biased but much less so...)

hist(comp$h2)
mean(comp$h2)
mean(comp$h2) - 2 * sqrt(var(comp$h2) / nrow(comp))
# for genetic variation of about 0.5 this seems off.

runner %>% select(contains('h2')) %>% plot()

ggplot(runner, aes(x = lande6, y = d2)) + 
  geom_point() +
  geom_segment(aes(x = 2, xend = 3, y = 2, yend = 3), colour = 'blue')
# all of these suuuuuck

ggplot(runner, aes(x = lande7, y = d2)) + 
  geom_point() +
  geom_segment(aes(x = 2, xend = 3, y = 2, yend = 3), colour = 'blue')

ggplot(runner, aes(x = lande8, y = d2)) + 
  geom_point() +
  geom_segment(aes(x = 2, xend = 3, y = 2, yend = 3), colour = 'blue')

# plot the landes
runner %>% select(c('d2', contains('lande'))) %>% plot()

### Conclusion of above: except for one that's obviously wrong (??),
# can't find a way to make the lande expression match up with sim outputs

# so, look further upstream
# look at zwt

lis2 = vector('list', 1000)

set.seed(1572)

for (i in 1:1000) {
  sim.output = sim( a = c(1/2, -1/2), params = pars[i,])
  
  demo.summ = sim.output %>%
    group_by(gen) %>%
    summarise(n = n(),
              gbar = mean(g_i),
              gvar = var(g_i),
              zbar = mean(z_i),
              zvar = var(z_i),
              evar = var(z_i - g_i),
              zwt  = mean(z_i * w_i) / mean(w_i))
  
  gene.summ = sim.output %>%
    select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
    gather(key = loc.copy, value = val, -gen) %>%
    mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
    group_by(gen, locus) %>%
    summarise(p = mean(val > 0)) %>%
    group_by(gen) %>%
    summarise(v = sum(2 * p * (1 - p)) / pars$n.loci[1])
  
  lis2[[i]] = cbind(demo.summ, gene.summ %>% select(-gen)) 
  
  print(paste0(i, ' of ', nrow(pars)))
}

pops2 = unroller(lis2)

com2 = pops2 %>%
  mutate(th = pars$theta[1]) %>%
  group_by(trial) %>%
  mutate(d1 = th - zbar[1],
         d2 = th - zbar[2],
         # fitness-weighted phenotypes straight from simulations
         dw1a = th - zwt,
         # lande's equation, p. 322, between (13) and (14)
         dw1b = d1 * (pars$wfitn[1]^2) / (pars$wfitn[1]^2 + zvar)) %>%
  ungroup() %>%
  filter(gen < 2)

# Differences between dw1 and dw2...?
com2 %>%
  ggplot() +
  geom_point(aes(x = dw1a, y = dw1b)) +
  geom_segment(aes(x = 2, xend = 2.5, y = 2, yend = 2.5), colour = 'blue')
# Yes, lande's expression is correct (should just be the integral)

com2 %>%
  mutate(wsum = (d1 * evar + dw1a * gvar) / (evar + gvar)) %>%
  ggplot() +
  geom_point(aes(x = wsum, y = d2)) +
  geom_segment(aes(x = 2, xend = 3, y = 2, yend = 3), colour = 'blue')

# yep - this weighted sum looks just like the lande eq. for zwt

com2 %>%
  mutate(wsum = (d1 * evar + dw1a * gvar) / (evar + gvar)) %>%
  ggplot() +
  geom_point(aes(x = wsum, y = d2)) +
  geom_segment(aes(x = 2, xend = 3, y = 2, yend = 3), colour = 'blue')

com2 %>%
  mutate(h2 = gvar / (gvar + evar),
         wzum = dw1a * h2 + d1 * (1-h2)) %>%
  ggplot() +
  geom_point(aes(x = wzum, y = d2)) +
  geom_segment(aes(x = 2, xend = 2.75, y = 2, yend = 2.75), colour = 'blue')

# In both cases, seeing the same pattern
# the new phenotype is *greater* than it the heritability-weighted sum
# (this means that either the expression is wrong, or h2 is wrong).

# So what's the problem...?
# Think mechanistically...

# the relized d2 is closer to theta, i.e., z2 is smaller than expected,
# i.e., k is *larger* than expected
# i.e., (w^2 + sig.e^2) / (w^2 + sig.e^2 + sig.a^2) is *larger* than expected
# i.e., (1 - h2 (sig^2) / (w^2 + sig^2)) is *larger* than expected
# (could be: h2, sig.e, sig, sig.a)

### Go back to the big batch of sims.

# What is the correct k?

comp %>%
  mutate(k_obs = d2 / d1) %>%
  ggplot() +
  geom_histogram(aes(x = k_obs))

comp %>%
  mutate(w2 = pars$wfitn[1]^2,
         k_obs = d2 / d1) %>%
  ggplot() +
  geom_point(aes(x = (w2 + evar) / (w2 + zvar),
                 y = k_obs)) +
  geom_segment(aes(x = 0.8, y = 0.8, xend = 1.0, yend = 1.0),
               colour = 'blue') +
  stat_smooth(aes(x = (w2 + evar) / (w2 + zvar),
                  y = k_obs),
              colour = 'orange',
              method = 'lm')

comp %>%
  mutate(w2 = pars$wfitn[1]^2,
         k_obs = d2 / d1) %>%
  ggplot() +
  geom_point(aes(x = (w2 + evar) / (w2 + evar + gvar),
                 y = k_obs)) +
  geom_segment(aes(x = 0.8, y = 0.8, xend = 1.0, yend = 1.0),
               colour = 'blue') +
  stat_smooth(aes(x = (w2 + evar) / (w2 + zvar),
                  y = k_obs),
              colour = 'orange',
              method = 'lm')

# Here too, estimates are biased
# k_obs is greater than what is k as predicted above
# (because we just dropped the d1 term from both sides of the equation)

comp %>%
  mutate(w2 = pars$wfitn[1]^2,
         k_obs = d2 / d1) %>%
  ggplot() +
  geom_point(aes(x = 1 - (evar / zvar) * (zvar / (w2 + zvar)),
                 y = k_obs)) +
  geom_segment(aes(x = 0.8, y = 0.8, xend = 1.0, yend = 1.0),
               colour = 'blue') +
  stat_smooth(aes(x = (w2 + evar) / (w2 + zvar),
                  y = k_obs),
              colour = 'orange',
              method = 'lm')

comp %>% ggplot() + geom_point(aes(x = v, y = gvar), alpha = 0.2)
mean(comp$gvar)
# these are the same on average but why aren't they 1-1...?
# and why is thre such hideously massive variation in gvar?

# (although neither of these seem to be causing the problem)

### Heritability - maybe the model has heritability wrong

# Load in pedigreed sims
# (note: overrides sim functions)
source('base_model_source/sim_functions_pedigreed.R')

set.seed(409)
herit.sim = sim.par.id(params = pars[1,])

head(herit.sim)

midp = merge(x = herit.sim %>% filter(gen %in% 2),
      y = herit.sim %>% filter(gen %in% 1) %>% select(i, g_i, z_i),
      # x is offspring, y is parent
      by.x = 'mom_i', by.y = 'i', suffixes = c('', '_mom')) %>%
  merge(y = herit.sim %>% filter(gen %in% 1) %>% select(i, g_i, z_i),
        by.x = 'dad_i', by.y = 'i', suffixes = c('', '_dad')) %>%
  mutate(z_midp = (z_i_mom + z_i_dad) / 2,
         g_midp = (g_i_mom + g_i_dad) / 2)

midp

with(midp, cov(z_midp, z_i) / var(z_i)) # 0.229?? that's very small
# but is this the same thing as broad sense?
with(midp, var(g_i) / var(z_i)) # 0.598...
with(midp, cov(g_i, z_i) / var(z_i)) # 0.579...
# this is just one trial.

# Run a loop.

herit.trials = 2
herit.list = vector('list', length = herit.trials)

set.seed(578)
for (i in 1:herit.trials) herit.list[[i]] = sim.par.id(params = pars[i,])

herit.sims = unroller(herit.list)

with(herit.sims, table(trial, gen))

# output should be nrow = 49 + 47

# Can I just include trial in merge...?
all.midps = merge(x = herit.sims %>% filter(gen %in% 2),
      y = herit.sims %>% filter(gen %in% 1) %>% select(i, g_i, z_i, trial),
      # x is offspring, y is parent
      by.x = c('mom_i', 'trial'), by.y = c('i', 'trial'), suffixes = c('', '_mom')) %>%
  merge(y = herit.sims %>% filter(gen %in% 1) %>% select(i, g_i, z_i, trial),
        by.x = c('dad_i', 'trial'), by.y = c('i', 'trial'), suffixes = c('', '_dad')) %>%
  mutate(z_midp = (z_i_mom + z_i_dad) / 2,
         g_midp = (g_i_mom + g_i_dad) / 2)

table(all.midps$trial)
# hell yes.

# Try a bigger batch of sims

herit.trials = 500
herit.list = vector('list', length = herit.trials)

set.seed(578)
for (i in 1:herit.trials) herit.list[[i]] = sim.par.id(params = pars[i,])

herit.sims = unroller(herit.list)

all.midps = merge(x = herit.sims %>% filter(gen %in% 2),
                  y = herit.sims %>% filter(gen %in% 1) %>% select(i, g_i, z_i, trial),
                  # x is offspring, y is parent
                  by.x = c('mom_i', 'trial'), by.y = c('i', 'trial'), suffixes = c('', '_mom')) %>%
  merge(y = herit.sims %>% filter(gen %in% 1) %>% select(i, g_i, z_i, trial),
        by.x = c('dad_i', 'trial'), by.y = c('i', 'trial'), suffixes = c('', '_dad')) %>%
  mutate(z_midp = (z_i_mom + z_i_dad) / 2,
         g_midp = (g_i_mom + g_i_dad) / 2)

table(all.midps$trial)

all.h2s = all.midps %>%
  group_by(trial) %>%
  mutate(var.a = var(g_i_dad) + var(g_i_mom),
         cov.mpz = cov(z_midp, z_i),
         narrow.h2.var = var.a / var(z_i),
         broads.h2.var = var(g_i) / var(z_i),
         broads.h2.cov = (cov(g_i, z_i) / (sd(g_i) * sd(z_i)))^2,
         broads.h2.cor = cor(g_i, z_i)^2,
         cov.g.z = cov(g_i, z_i),
         var.g = var(g_i)) %>%
  ungroup()

all.h2s %>% select(contains('h2')) %>% plot()
all.h2s %>% select(contains('h2')) %>% cor()
# das heck? correlation here is quite bad.
# (except that squared correlation IS the cov. expression)
# so - the var expression is _not_ the squared correlation
# interesting.

# is var(g) same as cov(g, z)?
plot(cov.g.z ~ var.g, all.h2s)
abline(0, 1, col = 'blue')
# plausibly is the same?
lm(cov.g.z ~ var.g, all.h2s) %>% summary()
abline(lm(cov.g.z ~ var.g, all.h2s), col = 'orange')
# This actually looks good
# some evidence for non-zero intercept, non-1 slope
# hmm...
lm(cov.g.z ~ var.g - 1, all.h2s) %>% summary()
abline(lm(cov.g.z ~ var.g - 1, all.h2s), col = 'red')
# still slightly less than 1.
# This looks like there is the tiniest amount of negative covariation between g and e??

with(all.h2s, hist(cov.g.z / var.g))
with(all.h2s, mean(cov.g.z / var.g))
# covariance is *slightly* smaller

all.h2s %>% select(contains('h2')) %>% apply(2, mean)
# broad sense heritabilities look fine
# very weird that the covariance-based one is lower (also that correlation is basically zero!)
# narrow sense heritabilities are g-darn whack though
