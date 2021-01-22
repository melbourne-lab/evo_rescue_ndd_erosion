library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

source('base_model_source/sim_functions.R')

pars = data.frame(
  end.time = 15,
  init.row = 500,
  n.pop0 = 30,
  n.loci = 25,
  w.max = 2,
  theta = 2.75,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0
)

set.seed(992)

sim1 = sim(params = pars)

sim1 %>% group_by(gen) %>% summarise(n = n(), wbar = mean(w_i))

new.gen = sim1 %>%
  filter(gen %in% max(gen)) %>%
  sample_n(10, replace = TRUE) %>%
  mutate(fem = sample(rep(c(TRUE, FALSE), 5)),
         r_i = rpois(n = 10, lambda = ifelse(fem, 2 * w_i * exp(-pars$alpha * 10), 0)))

propagate.sim(params = pars, popn = new.gen)
         

# Actually, what if we just removed the stochasticity? Resample sex and resample r_i?

new.gen = sim1 %>%
  filter(gen %in% max(gen)) %>%
  mutate(fem = sample(size = nrow(.), x = c(TRUE, FALSE), replace = TRUE),
         r_i = rpois(n = nrow(.), lambda = ifelse(fem, 2 * w_i * exp(-pars$alpha * nrow(.)), 0)))

while (!sum(new.gen$r_i)) {
  new.gen = new.gen %>%
    mutate(fem = sample(size = nrow(.), x = c(TRUE, FALSE), replace = TRUE),
           r_i = rpois(n = nrow(.), lambda = ifelse(fem, 2 * w_i * exp(-pars$alpha * nrow(.)), 0)))
}

new.gen

propagate.sim(params = pars, popn = new.gen)

### Okay, define some functions here

revive.popn = function(popn, params) {
  # Function to revive populations
  
  alph = params$alpha
  
  revived.popn = popn %>%
    filter(gen %in% max(gen)) %>%
    mutate(fem = sample(size = nrow(.), x = c(TRUE, FALSE), replace = TRUE),
           r_i = rpois(n = nrow(.), lambda = ifelse(fem, 2 * w_i * exp(-alph * nrow(.)), 0)))
  
  while(!sum(revived.popn$r_i)) {
    revived.popn = revived.popn %>%
      filter(gen %in% max(gen)) %>%
      mutate(fem = sample(size = nrow(.), x = c(TRUE, FALSE), replace = TRUE),
             r_i = rpois(n = nrow(.), lambda = ifelse(fem, 2 * w_i * exp(-alph * nrow(.)), 0)))
  }
  
  return(revived.popn)
  
}

# Functon to revive sims.

finish.revive.sim = function(a, popn, params) {
  
  while(max(popn$gen) < params$end.time) {
    
    
    
  }
  
}