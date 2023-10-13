# Robustness check for environmental contribution to phenotypic variance
# (sigma^2_e in model)
# SN - run 8 Sept. 2023

library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

rm(list = ls())

source('sim_functions.R')

trials.per = 500

pars = expand.grid(
  alpha = c(0, 0.0035),
  n.pop0 = c(20, 100),
  lowvar = c(TRUE, FALSE),
  sig.e = sqrt(c(0.25, 1, 2)),
  trial = 1:trials.per
) %>%
  mutate(
    size.thresh = 50000,
    init.row = 1e6,
    end.time = 16,
    n.loci = 25,
    w.max = 2,
    wfitn = sqrt(3.5),
    theta = 2.8,
    pos.p = 0.5
  )

set.seed(4290) 

sim.out = mclapply(
  pars %>% mutate(try.no = 1:nrow(.)) %>% split(.$try.no),
  function (params) {
    if (params$lowvar) {
      sim.raw = sim(
        params, 
        init.popn = init.sim(params) %>%
          mutate_at(paste0('a', 1:6), function(x) -1/2) %>%
          mutate_at(paste0('b', 1:6), function(x) -1/2) %>%
          mutate_at(paste0('a', 7:12), function(x) 1/2) %>%
          mutate_at(paste0('b', 7:12), function(x) 1/2)
      )
    } else {
      sim.raw = sim(params)
    }
    cat(params$trial / trials.per)
    cat(" ")
    sim.raw %>%
      group_by(gen) %>%
      summarise(n = n()) %>%
      mutate(
        alpha = params$alpha,
        n.pop0 = params$n.pop0,
        lowvar = params$lowvar,
        sig.e  = params$sig.e,
        trial = params$try.no
      )
  },
  mc.cores = 16
) %>%
  do.call(what = rbind) %>%
  mutate(gen = gen - 1)

write.csv(sim.out,'simulations/outputs/robust_checks_out/sig2e_robust_check.csv', row.names = FALSE)
