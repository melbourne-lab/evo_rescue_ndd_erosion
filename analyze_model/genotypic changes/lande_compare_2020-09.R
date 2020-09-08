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
  mutate(wsum = (d1 * evar + dw1a * gvar) / (evar + gvar))) %>%
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

# So what's the problem...?
# Think mechanistically...

# the relized d2 is closer to theta, i.e., z2 is smaller than expected,
# i.e., k is *larger* than expected
# i.e., (w^2 + sig.e^2) / (w^2 + sig.e^2 + sig.a^2) is *larger* than expected
# i.e., (1 - h2 (sig^2) / (w^2 + sig^2)) is *larger* than expected
# (could be: h2, sig.e, sig, sig.a)
