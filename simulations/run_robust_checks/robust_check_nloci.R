# Robustness check on number of loci
# SN run 12 Sept 2023

library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

rm(list = ls())

source('sim_functions.R')

trials.per = 1000

pars = expand.grid(
  alpha = c(0, 0.0035),
  n.pop0 = c(20, 100),
  lowvar = c(TRUE, FALSE),
  n.loci = c(9, 36, 100),
  trial = 1:trials.per
) %>%
  mutate(
    size.thresh = 50000,
    init.row = 1e6,
    end.time = 16,
    w.max = 2,
    wfitn = sqrt(3.5),
    theta = 2.8,
    sig.e = sqrt(0.5),
    pos.p = 0.5
  )

set.seed(900491) 

sim.out.m100 = mclapply(
  pars %>% filter(n.loci %in% 100) %>% mutate(try.no = 1:nrow(.)) %>% split(.$try.no),
  function (params) {
    if (params$lowvar) {
      sim.raw = sim(
        params, 
        init.popn = init.sim(params) %>%
          mutate_at(paste0('a', 1:25), function(x) -1/2) %>%
          mutate_at(paste0('b', 1:25), function(x) -1/2) %>%
          mutate_at(paste0('a', 26:50), function(x) 1/2) %>%
          mutate_at(paste0('b', 26:50), function(x) 1/2)
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
        n.loci = params$n.loci,
        trial = params$try.no
      )
  },
  mc.cores = 20
) %>%
  do.call(what = rbind)

sim.out.m36 = mclapply(
  pars %>% filter(n.loci %in% 36) %>% mutate(try.no = 1:nrow(.)) %>% split(.$try.no),
  function (params) {
    if (params$lowvar) {
      sim.raw = sim(
        params, 
        init.popn = init.sim(params) %>%
          mutate_at(paste0('a', 1:9), function(x) -1/2) %>%
          mutate_at(paste0('b', 1:9), function(x) -1/2) %>%
          mutate_at(paste0('a', 10:18), function(x) 1/2) %>%
          mutate_at(paste0('b', 10:18), function(x) 1/2)
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
        n.loci = params$n.loci,
        trial = params$try.no
      )
  },
  mc.cores = 20
) %>%
  do.call(what = rbind)

sim.out.m9 = mclapply(
  pars %>% filter(n.loci %in% 9) %>% mutate(try.no = 1:nrow(.)) %>% split(.$try.no),
  function (params) {
    if (params$lowvar) {
      sim.raw = sim(
        params, 
        init.popn = init.sim(params) %>%
          mutate_at(paste0('a', 1:2), function(x) -1/2) %>%
          mutate_at(paste0('b', 1:2), function(x) -1/2) %>%
          mutate_at(paste0('a', 3:4), function(x) 1/2) %>%
          mutate_at(paste0('b', 3:4), function(x) 1/2)
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
        n.loci = params$n.loci,
        trial = params$try.no
      )
  },
  mc.cores = 20
) %>%
  do.call(what = rbind)

write.csv(
  rbind(sim.out.m100, sim.out.m36, sim.out.m9) %>% mutate(gen = gen - 1),
  'simulations/outputs/robust_checks_out/nloci_robust_check.csv', 
  row.names = FALSE
)
