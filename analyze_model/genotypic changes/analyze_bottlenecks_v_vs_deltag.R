# Initial run throughs of main simulations suggest little effect of
# bottlenecking on population dynamics. What gives?
# Here I run trials (evolving, alpha = 0) and look at relationships between the
# genetic variance, change in genotype, and population growth rate changes.
# SN - July 14 2020

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

# Load source populations
# (do this now to get the population numbers)
neck.popns = read.csv('simulations/prepare_simulations/outputs/all_k_pops_min_gen.csv',
                      stringsAsFactors = FALSE)

# Define number of trials per source populations
trials.per = 26

pars = expand.grid(pop.id = unique(neck.popns$trial),
                   pop.trial = 1:trials.per,
                   n.pop0 = c(20, 100)) %>%
  mutate(global.trial = 1:nrow(.),
         end.time = 2,
         init.row = 1e4,
         n.loci = 25,
         w.max = 2,
         theta = 2.6,
         wfitn = sqrt(1 / 0.14 / 2),
         sig.e = 0,
         pos.p = 0.5,
         alpha = 0)

liszt.vlo = vector('list', nrow(pars))

set.seed(7756758)

for (i in 1:nrow(pars)) {
  popid = pars$pop.id[i]
  n0 = pars$n.pop0[i]
  liszt.vlo[[i]] = sim( a = c(1/2, -1/2), 
                        params = pars[i,],
                        init.popn = neck.popns %>%
                          filter(trial %in% popid ) %>% # [neck.popns$trial %in% pars$pop.id[i],] %>%
                          sample_n(n0, replace = FALSE) %>%
                          select(-trial)
  )
  print(paste0('bottle ', i, ' of ', nrow(pars)))
}

all.vlo.trials = unroller(liszt.vlo)

# head(all.vlo.trials)
# 
# all.vlo.trials %>%
#   group_by(trial, gen) %>%
#   summarise(n = n())

all.vlo.n = all.vlo.trials %>%
  group_by(trial, gen) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  rbind(expand.grid(trial = 1:nrow(pars), gen = 1:pars$end.time[1], n = 0)) %>%
  group_by(trial, gen) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  merge(y = pars, by.x = 'trial', by.y = 'global.trial')

### Normie populations (i.e., non-bottlenecked)

liszt.vhi = vector('list', nrow(pars))

set.seed(69594)

for (i in 1:nrow(pars)) {
  liszt.vhi[[i]] = sim( params = pars[i,])
  print(paste0('non-bottle ', i, ' of ', nrow(pars)))
}

all.vhi.trials = unroller(liszt.vhi)

# head(all.vlo.trials)
# 
# all.vhi.trials %>%
#   group_by(trial, gen) %>%
#   summarise(n = n())

all.vhi.n = all.vhi.trials %>%
  group_by(trial, gen) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  rbind(expand.grid(trial = 1:nrow(pars), gen = 1:pars$end.time[1], n = 0)) %>%
  group_by(trial, gen) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  merge(y = pars, by.x = 'trial', by.y = 'global.trial')

### Combine the results.

all.sizes = 
  rbind(
    all.vlo.n %>% mutate(bottleneck = TRUE),
    all.vhi.n %>% mutate(bottleneck = FALSE, trial = trial + max(all.vlo.n$trial))
  )

all.demos = rbind(
  all.vlo.trials %>% 
    mutate(bottleneck = TRUE) %>%
    merge(y = pars %>% select(global.trial, n.pop0),
          by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = TRUE),
  all.vhi.trials %>%
    merge(y = pars %>% select(global.trial, n.pop0),
          by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = FALSE, trial = trial + max(all.vlo.n$trial))
) %>%
  group_by(trial, gen, bottleneck, n.pop0) %>%
  summarise(gbar = mean(g_i),
            gvar = var(g_i),
            zbar = mean(z_i),
            zvar = var(z_i),
            wbar = mean(w_i),
            wvar = var(w_i))

all.genes = rbind(
  all.vlo.trials %>% 
    mutate(bottleneck = TRUE) %>%
    merge(y = pars %>% select(global.trial, n.pop0),
          by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = TRUE),
  all.vhi.trials %>%
    merge(y = pars %>% select(global.trial, n.pop0),
          by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = FALSE, trial = trial + max(all.vlo.n$trial))
) %>%
  select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(gen, trial, bottleneck, n.pop0)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(trial, bottleneck, n.pop0, gen, locus) %>%
  summarise(p = mean(val > 0)) %>%
  group_by(trial, bottleneck, n.pop0, gen) %>%
  summarise(m.fix.pos = sum(p == 1),
            m.fix.neg = sum(p == 0),
            v = sum(2 * p * (1 - p)) / pars$n.loci[1]) #%>%
  #group_by(bottleneck, n.pop0, gen) %>%
  summarise(p.fix.pos = mean(p.fix.pos),
            p.fix.neg = mean(p.fix.neg),
            vbar = mean(v),
            vvar = var(v))

all.delts = 
  merge(
  x = all.genes %>%
    group_by(trial, bottleneck, n.pop0) %>%
    summarise(del.v = diff(v),
              v1 = v[1]) %>%
    ungroup(),
  y = all.demos %>%
    group_by(trial, bottleneck, n.pop0) %>%
    summarise(del.g = diff(gbar),
              g1 = gbar[1]) %>%
    ungroup(),
  by = c('trial', 'bottleneck', 'n.pop0')
) %>%
  merge(y = all.sizes %>%
          group_by(trial, bottleneck, n.pop0) %>%
          summarise(log.l = log(n[2] / n[1])) %>%
          ungroup(),
        by = c('trial', 'bottleneck', 'n.pop0'))

ggdelt = ggplot(all.delts)

# Question 1: does v_1 predict delta_g?
ggdelt +
  geom_point(aes(x = v1, y = del.g / (g1 - 2.6),
                 colour = bottleneck,
                 shape = factor(n.pop0))) +
  stat_smooth(aes(x = v1, y = del.g / (g1 - 2.6)),
              method = 'lm')
# this is actually plausible... 
# (although variance would depend on n.pop0)

ggdelt +
  geom_point(aes(x = v1, y = del.g / (g1 - 2.6))) +
  stat_smooth(aes(x = v1, y = del.g / (g1 - 2.6)),
              method = 'lm') +
  facet_wrap(bottleneck ~ n.pop0,
             scales = "free_x")

# The relationship here is weak... but I think it exists

all.delts %>% 
  mutate(prop.dg = del.g / (g1-2.6)) %>% 
  lm(formula = prop.dg ~ v1) %>% 
  summary()

all.delts %>% 
  mutate(prop.dg = del.g / (g1-2.6)) %>% 
  lm(formula = prop.dg ~ v1 - 1) %>% 
  summary()

# variance needs to be adjusted for n.pop0
# but this relationship is strong
# (I forget - is R2 interpretable for zero-int regression?
# that thing is massive... probably total baloney)
# R2 is low for the other one

ggdelt +
  geom_point(aes(x = v1, y = del.g)) +
  stat_smooth(aes(x = v1, y = del.g),
              method = 'lm') +
  facet_wrap(bottleneck ~ n.pop0,
             scales = "free_x")
# this is bad

# Question 2: does V predict log.l?
ggdelt +
  geom_point(aes(x = v1, y = log.l)) +
  stat_smooth(aes(x = v1, y = log.l),
              method = 'lm') +
  facet_wrap(bottleneck ~ n.pop0,
             scales = "free_x")
# the relationship here is... not very strong

ggdelt +
  geom_point(aes(x = v1, y = log.l,
                 colour = bottleneck)) +
  stat_smooth(aes(x = v1, y = log.l),
              method = 'lm') +
  facet_wrap(~n.pop0)
# this looks very bad

# Question 3: does delta_g correlate at all with log.l?
ggdelt +
  geom_point(aes(x = del.g / (g1 - 2.6), y = log.l)) +
  stat_smooth(aes(x = del.g / (g1 - 2.6), y = log.l),
              method = 'lm') +
  facet_wrap(bottleneck ~ n.pop0,
             scales = "free_x")
# also no
# hmm...

ggdelt +
  geom_point(aes(x = del.g / (g1 - 2.6), y = log.l,
                 colour = bottleneck)) +
  stat_smooth(aes(x = del.g / (g1 - 2.6), y = log.l),
              method = 'lm') +
  facet_wrap(~ n.pop0)
# residuals here are definitely a problem

all.delts %>% 
  mutate(prop.dg = del.g / (g1-2.6)) %>% 
  filter(n.pop0 > 50) %>%
  lm(formula = log.l ~ prop.dg) %>% 
  summary()
# p-vals are non-significant but of course resids are bad
# R2 is also basically zero

# Looking at the distribution of genotypic changes

ggdelt + 
  geom_histogram(aes(x = del.g / (g1 - 2.6)), binwidth = 0.025) +
  facet_wrap(bottleneck ~ n.pop0)

ggdelt + 
  geom_density(aes(x = del.g, 
                   group = interaction(n.pop0, bottleneck),
                   colour = interaction(n.pop0, bottleneck))) +
  scale_colour_manual(values = c('red', 'pink', 'blue', 'lightblue'))

ggdelt + 
  geom_density(aes(x = del.g / (g1 - 2.6), 
                   group = interaction(n.pop0, bottleneck),
                   colour = interaction(n.pop0, bottleneck))) +
  scale_colour_manual(values = c('red', 'pink', 'blue', 'lightblue'))
# does seem like you get _slightly_ less change in g with bottleneck
# (although at n0 = 20 it's pretty tough to see - lots of varance)

ggdelt + 
  geom_density(aes(x = log.l, 
                   group = interaction(n.pop0, bottleneck),
                   colour = interaction(n.pop0, bottleneck))) +
  scale_colour_manual(values = c('red', 'pink', 'blue', 'lightblue'))
# but there isn't much of an effect on log lambda at this point

ggdelt +
  geom_density(aes(x = del.v, 
                   group = interaction(n.pop0, bottleneck),
                   colour = interaction(n.pop0, bottleneck))) +
  scale_colour_manual(values = c('red', 'pink', 'blue', 'lightblue'))
# so there's slightly more loss of variation in non-bottled populations

# but if we adjust by initial amount of variation,
ggdelt +
  geom_density(aes(x = del.v/v1, 
                   group = interaction(n.pop0, bottleneck),
                   colour = interaction(n.pop0, bottleneck))) +
  scale_colour_manual(values = c('red', 'pink', 'blue', 'lightblue'))
# differences look like they disappear
# so in a sense the *proportional* loss of variation is the same

# What does this mean??
# - There is a relationship between v1 and delta.g
# Decreasing the variance (v) decreases the genotypic change
# (the genotypes approach the optimum at a slower pace)
# This is in accordance with the Lande model
# (not sure yet if the approach is geometric though)
# - The relationship between v1 log_lambda is weak
# (presumably because there is extra phenotypic noise obscuring)
# meaning that genetic variance has close to no impact on actual
# population growth
# this explains the decoupling
# - Likewise there is little relationship between g and log-lamb
# (probably also due to phenotpyic noise)

# Try this next with no phenotypic noise

pars2 = pars %>% mutate(sig.e = 0)

list2.vlo = vector('list', nrow(pars2))

set.seed(7756758)

for (i in 1:nrow(pars2)) {
  popid = pars2$pop.id[i]
  n0 = pars2$n.pop0[i]
  list2.vlo[[i]] = sim( a = c(1/2, -1/2), 
                        params = pars2[i,],
                        init.popn = neck.popns %>%
                          filter(trial %in% popid ) %>% # [neck.popns$trial %in% pars$pop.id[i],] %>%
                          sample_n(n0, replace = FALSE) %>%
                          select(-trial)
  )
  print(paste0('bottle ', i, ' of ', nrow(pars)))
}

vlo.trials.g = unroller(list2.vlo)

### non-bottlenecked populations

list2.vhi = vector('list', nrow(pars2))

set.seed(69594)

for (i in 1:nrow(pars)) {
  list2.vhi[[i]] = sim( params = pars2[i,])
  print(paste0('non-bottle ', i, ' of ', nrow(pars)))
}

vhi.trials.g = unroller(list2.vhi)

# Summary stats

all.demos.g = rbind(
  vlo.trials.g %>% 
    mutate(bottleneck = TRUE) %>%
    merge(y = pars %>% select(global.trial, n.pop0),
          by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = TRUE),
  vhi.trials.g %>%
    merge(y = pars %>% select(global.trial, n.pop0),
          by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = FALSE, trial = trial + max(all.vlo.n$trial))
) %>%
  group_by(trial, gen, bottleneck, n.pop0) %>%
  summarise(gbar = mean(g_i),
            gvar = var(g_i),
            zbar = mean(z_i),
            zvar = var(z_i),
            wbar = mean(w_i),
            wvar = var(w_i))

all.genes.g = rbind(
  vlo.trials.g %>% 
    mutate(bottleneck = TRUE) %>%
    merge(y = pars %>% select(global.trial, n.pop0),
          by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = TRUE),
  vhi.trials.g %>%
    merge(y = pars %>% select(global.trial, n.pop0),
          by.x = 'trial', by.y = 'global.trial') %>%
    mutate(bottleneck = FALSE, trial = trial + max(all.vlo.n$trial))
) %>%
  select(-c(i, g_i, z_i, w_i, r_i, fem)) %>%
  gather(key = loc.copy, value = val, -c(gen, trial, bottleneck, n.pop0)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(trial, bottleneck, n.pop0, gen, locus) %>%
  summarise(p = mean(val > 0)) %>%
  group_by(trial, bottleneck, n.pop0, gen) %>%
  summarise(m.fix.pos = sum(p == 1),
            m.fix.neg = sum(p == 0),
            v = sum(2 * p * (1 - p)) / pars$n.loci[1])

all.delts.g = 
  merge(
    x = all.genes.g %>%
      group_by(trial, bottleneck, n.pop0) %>%
      summarise(del.v = diff(v),
                v1 = v[1]) %>%
      ungroup(),
    y = all.demos.g %>%
      group_by(trial, bottleneck, n.pop0) %>%
      summarise(del.g = diff(gbar),
                g1 = gbar[1]) %>%
      ungroup(),
    by = c('trial', 'bottleneck', 'n.pop0')
  )

ggdelt.g = ggplot(all.delts.g)

ggdelt.g +
  geom_point(aes(x = v1, y = del.g / (g1 - 2.6))) +
  stat_smooth(aes(x = v1, y = del.g / (g1 - 2.6)),
              method = 'lm') +
  facet_wrap(bottleneck ~ n.pop0,
             scales = "free_x")
# Here the relationship is also very very weak.

all.delts.g %>% 
  mutate(prop.dg = del.g / (g1-2.6)) %>% 
  lm(formula = prop.dg ~ v1) %>% 
  summary()

# Here the variance is still very large
# I.e., having no phenotypic variance doesn't have a huge effect.