# Robustness check on number of loci
# SN run 12 Sept 2023

library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

rm(list = ls())

source('sim_bh_functions.R')

trials.per = 1000

pars = expand.grid(
  alpha = c(0, 0.0035),
  n.pop0 = c(20, 100),
  lowvar = c(TRUE, FALSE),
  trial = 1:trials.per,
  delta = c(.5, .75, 1)
) %>%
  mutate(
    size.thresh = 10000,
    init.row = 1e6,
    end.time = 16,
    w.max = 2,
    n.loci = 25,
    wfitn = sqrt(3.5),
    theta = 2.8,
    sig.e = sqrt(0.5),
    pos.p = 0.5,
    wbar  = w.max * sqrt(wfitn^2 / (wfitn^2 + sig.e^2 + .5 - as.numeric(lowvar)*.25)),
    nstar =  log(wbar) / alpha,
    beta = ifelse(alpha > 0, (wbar^(1/delta) - 1) / nstar, 0)) %>%
  filter(!(!alpha & delta < 1))

set.seed(900491) 

sim.out = mclapply(
  pars %>% mutate(try.no = 1:nrow(.)) %>% split(.$try.no),
  function (params) {
    if (params$lowvar) {
      sim.raw = sim.bh(
        params, 
        init.popn = init.sim(params) %>%
          mutate_at(paste0('a', 1:6), function(x) -1/2) %>%
          mutate_at(paste0('b', 1:6), function(x) -1/2) %>%
          mutate_at(paste0('a', 7:12), function(x) 1/2) %>%
          mutate_at(paste0('b', 7:12), function(x) 1/2)
      )
    } else {
      sim.raw = sim.bh(params)
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
        delta = params$delta,
        beta  = params$beta,
        trial = params$try.no
      )
  },
  mc.cores = 16
) %>%
  do.call(what = rbind) %>%
  mutate(gen = gen - 1)

write.csv(
  sim.out,
  'simulations/outputs/robust_checks_out/bholt_robust_check.csv', 
  row.names = FALSE
)
