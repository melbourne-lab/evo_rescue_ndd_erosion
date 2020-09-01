# I'd love to be able to get summary statistics based on being extinct/extant in that timestep. This is absurdly difficult to do in dplyr. Doing it a very ugly and hacky way here. 

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


# Define number of trials per source populations
trials = 800

pars = data.frame(
  global.trial = 1:trials,
  end.time = 12,
  init.row = 1e4,
  n.loci = 25,
  n.pop0 = 60,
  w.max = 2,
  theta = 2.75,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0.0035)

liszt.vlo = vector('list', nrow(pars))

set.seed(8222)

for (i in 1:nrow(pars)) {
  pop.init = init.sim(a = c(1/2, -1/2),
                      params = pars[i,]) %>%
    mutate_at(paste0('a', 1:6), function(x) -1/2) %>%
    mutate_at(paste0('b', 1:6), function(x) -1/2) %>%
    mutate_at(paste0('a', 7:12), function(x) 1/2) %>%
    mutate_at(paste0('b', 7:12), function(x) 1/2)
  liszt.vlo[[i]] = sim( a = c(1/2, -1/2), 
                        params = pars[i,],
                        init.popn = pop.init
  ) %>%
    # Get population size in each step
    group_by(gen) %>%
    summarise(n = n()) %>%
    # Add zeros if needed
    rbind(data.frame(gen = 1:pars$end.time[1], n = 0)) %>%
    group_by(gen) %>%
    summarise(n = sum(n)) %>%
    ungroup()
  print(paste0('ndd lo var ', i, ' of ', nrow(pars)))
}

all.vlo.trials = unroller(liszt.vlo) %>%
  group_by(trial) %>%
  mutate(n.pop0 = n[1],
         ext.gen = max(gen[n>0]),
         extinct = !ext.gen %in% pars$end.time[1]) %>%
  ungroup()


eat.it = all.vlo.trials %>%
  merge(y = expand.grid(trial = 1:max(.$trial), ext.gen = min(.$ext.gen):12),
        by = 'trial', suffixes = c('', '.fake')) %>%
  group_by(trial) %>%
  filter(ext.gen > ext.gen.fake) %>%
  mutate(ext.gen = ifelse(ext.gen < 10, paste0('0', ext.gen), as.character(ext.gen)),
         ext.gen = paste0('ext', ext.gen),
         ext.gen.fake = ifelse(ext.gen.fake < 10, paste0('0', ext.gen.fake), ext.gen.fake),
         ext.gen.fake = paste0('ext', ext.gen.fake)) %>%
  ungroup() %>%
  gather(key = is.true.gen, val = ext.gen, ext.gen, ext.gen.fake) %>%
  distinct(trial, gen, ext.gen, .keep_all = TRUE) %>%
  mutate(is.true.gen = !grepl('fake', is.true.gen)) %>%
  spread(key = ext.gen, val = is.true.gen) 

eat.it %>%
  summarise_at(.funs = mean(n), .vars = 6:ncol(.))


# is the above data structure what I want?
# want to be able to go through every combination of T/F for ext03, ext04, etc...
# and calculate the mean of n

# ah fuck it!
# this is fucking stupid

bb = expand.grid(gen = 1:11, extinct = c(TRUE, FALSE), nbar = NA, nvar = NA)

for (timest in unique(bb$gen)) {
  this.step = all.vlo.trials %>% filter(gen %in% timest)
  bb$nbar[with(bb, gen %in% timest & extinct)] = mean(this.step$n[this.step$ext.gen %in% timest])
  bb$nvar[with(bb, gen %in% timest & extinct)] = var(this.step$n[this.step$ext.gen %in% timest])
  bb$nbar[with(bb, gen %in% timest & !extinct)] = mean(this.step$n[this.step$ext.gen > timest])
  bb$nvar[with(bb, gen %in% timest & !extinct)] = var(this.step$n[this.step$ext.gen > timest])
}

bb %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = nbar, group = extinct))
