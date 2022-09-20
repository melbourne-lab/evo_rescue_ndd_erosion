### Functions needed for running simulations


library(ggplot2) # only used for examining test results
library(dplyr)   # needed - version 1.0.7 used
library(tidyr)   # needed for spread/gather - version 1.1.3 used
library(tidyselect)  # needed for column selecting (alleles) - version 1.1.1 used

### Auxiliary functions

# Used for easily setting column names (couldn't find an easy way of doing this
# in tidyverse)
set.names = function(df, name.array) {
  names(df) = name.array
  df
}

# Used for adding rows to data frame
# (using pre-allocation here for expedience)
dim.add = function(df, rows, addition) {
  if (sum(is.na(df$i)) < nrow(addition)) {
    df = df %>%
      rbind(df %>%
              sample_n(size = rows, replace = TRUE) %>%
              mutate_all(function(x) NA))
  }
  if (nrow(addition)) {
    df[which.max(is.na(df$i)) + 0:(nrow(addition)-1),] = addition
  }
  return(df)
}

# Function for turning lists into data frame with trial info added
# (in retrospect adding trial here is not the best way to do this
# - would have been better to add a trial column directly after the sim)
unroller = function(sim.list) {
  sim.list %>%
    do.call(what = rbind) %>%
    mutate(trial = cumsum(c(1, as.numeric(diff(gen) < 0))))
}

### Main simulation functions

# Wrapper for initializing a population for simulation
init.sim = function(a = c(1/2, -1/2), params) {
  
  # Inputs:
  # a - an array of length two (bi-allelic model)
  # each element is a contribution to the genotype
  # n.b. default is simplest case: -1/2, 1/2
  # params - a data frame of parameters
  
  n.loci = params$n.loci
  # number of loci determining the genotype
  n.pop0 = params$n.pop0
  # initial population size
  w.max = params$w.max
  # max population size
  theta = params$theta
  # optimal phenotypic value
  wfitn = params$wfitn
  # standard deviation of the selection pressure
  sig.e = params$sig.e
  # standard dev. of environmental phenoytpic noise
  p.pos = ifelse('pos.p' %in% names(params), params$pos.p, 0.5)
  # Initial frequency of the positive allele
  alpha = ifelse('alpha' %in% names(params), params$alpha, 0)
  # strength of density dependence
  # (set to zero if not provided)
    
  
  # A character (string) array for handy indexing
  names.array = paste0(c('a', 'b'), rep(1:n.loci, each = 2))
  
  # Initialize population.
  # matrix() - generate a matrix of alleles, with one row for each
  # individual, with two copies for each individual
  # sampling from a at an equal allele frequency
  # data.frame() - convert the matrix into a data frame
  # set.names() - set names for the data frame to denote
  # that each column corresponds to one allele
  # mutate() - for each individual, calculate the following:
  #   - i:  a unique identifier for each individual
  #   - fem: assign sex to each individual
  #   - g_i: genetic value; sum of alleles, scaled by sqrt(n)
  #   (scaling by sqrt(n) ensures sums don't explode to infinity)
  #   - z_i: phenotype; genotypic value plus normally-distributed noise
  #   - r_i: number of offspring for each individual
  #   Poisson distributed, with mean
  #   W_i = W_max * exp(-k * (z_i - theta)^2)
  init.popn = matrix(sample(a, prob = c(p.pos, 1 - p.pos),
                            size = 2 * n.loci * n.pop0, replace = TRUE), 
                     nrow = n.pop0, 
                     ncol = 2 * n.loci) %>%
    data.frame() %>%
    set.names(name.array = all_of(names.array)) %>%
    mutate(g_i = apply(., 1, sum) / sqrt(n.loci),
           i = 1:n.pop0,
           fem = sample(c(TRUE, FALSE), size = n.pop0, replace = TRUE),
           z_i = rnorm(n.pop0, mean = g_i, sd = sig.e),
           w_i = w.max * exp(-(z_i - theta)^2 / (2 * wfitn^2)),
           r_i = rpois(n = n.pop0, lambda = ifelse(fem, 2 * w_i * exp(-alpha * n.pop0), 0)),
           gen = 1) %>%
    select(i, g_i, z_i, w_i, r_i, fem, gen, all_of(names.array))
  
  return(init.popn)
  
}

##### Tests of the above:
# (note: not run since early testing - can not guarantee quality of results)
# 
# pars = data.frame(n.loci = 20, n.pop0 = 40,
#                   w.max = 1.2, theta = 2.6,
#                   wfitn = sqrt(1 / 0.14),
#                   sig.e = 0.5)
# 
# popn0 = init.sim(a = c(-1/2, 1/2), params = pars)
# popn0
# 
# # Looks good.

# # Trying a sim with a custom allele frequency:
# 
# pars = data.frame(n.loci = 20, n.pop0 = 40,
#                   w.max = 1.2, theta = 2.6,
#                   wfitn = sqrt(1 / 0.14),
#                   sig.e = 0.5, pos.p = 0.5)
# 
# popn0 = init.sim(params = pars)
# popn0
# 
# pars = data.frame(n.loci = 20, n.pop0 = 40,
#                   w.max = 1.2, theta = 2.6,
#                   wfitn = sqrt(1 / 0.14),
#                   sig.e = 0.5, pos.p = 0.75)
# 
# popn0 = init.sim(params = pars)
# popn0
# popn0[, grep('^[ab]\\d', names(popn0), value = TRUE)] %>%
#   apply(2, function(x) mean(x > 0))
#
# # Try one with density dependence
# pars = data.frame(n.loci = 20, n.pop0 = 40,
#                   w.max = 1.2, theta = 2.6,
#                   wfitn = sqrt(1 / 0.14),
#                   sig.e = 0.5, alpha = 0.0035)
# 
# popna = init.sim(a = c(-1/2, 1/2), params = pars)
# popna
# hist(popn0$r_i)
# hist(popna$r_i)

# Wrapper for iterating the simulation forward oe timestep/generation
propagate.sim = function(a = c(1/2, -1/2), params, popn, evolve = TRUE) {
  
  # n.b. defaults are (1) alleles of same magnitude opposing effect and
  # (2) population should be evolving in each time step
  # (non-evolving treatment initializes population of same size as population
  # size in given timestep)
  
  n.loci = params$n.loci
  # number of loci determining the genotype
  n.pop0 = params$n.pop0
  # initial population size
  w.max = params$w.max
  # max population size
  theta = params$theta
  # optimal phenotypic value
  wfitn = params$wfitn
  # standard deviation of the selection pressure
  sig.e = params$sig.e
  # standard dev. of environmental phenoytpic noise
  alpha = ifelse('alpha' %in% names(params), params$alpha, 0)
  # strength of density dependence
  # (set to zero if not provided)
  
  # A character (string) array for handy indexing
  names.array = paste0(c('a', 'b'), rep(1:n.loci, each = 2))
  
  # a flag handling ensuring all the following are true:
  #   - any females? (else, can't have birth)
  #   - any males? (else, can't mate)
  #   - any offspring?
  if (any(popn$fem) & any(!popn$fem) & sum(popn$r_i[popn$fem])) {
    
    if (evolve) {
      
      offspring = cbind(
        # Maternal data frame:
        # Take the female rows in the data frame
        # Remove unnecessary columns (don't need to be inherited)
        # Rename columns to indicate alleles inhereted from mom
        #   NOTE: r_i also included here because we will need it later
        popn %>% 
          filter(fem) %>% 
          select(-c(i, g_i, w_i, z_i, fem, gen)) %>%
          set.names(paste(ifelse(grepl('^[ab]\\d', names(.)), 'mom', ''),
                          names(.), 
                          sep = '_')),
        # Paternal data frame
        # Sample these to get mating pairs, i.e.,
        #   draw from the pool of males once for each female
        # Select the paternal alleles
        # Rename columns to indicate alleles inhereited from dad
        # NOTE: this assumes that each mom mates with only one dad
        popn %>% 
          sample_n(size = sum(fem), 
                   weight = as.numeric(!fem) / sum(as.numeric(!fem)),
                   replace = TRUE) %>%
          select(all_of(names.array)) %>%
          set.names(paste('dad', names(.), sep = '_'))
      ) %>%
        # For each mating pair, duplicate by the number of offspring
        #   as determined by r_i
        uncount(weight = `_r_i`) %>%
        # Add a new column for unique ID of each individual
        #   (note - we'll have to remove this later for silly reasons)
        mutate(i = max(popn$i) + 1:nrow(.)) %>%
        # Use some cleverness to segregate alleles:
        #   create a row for each allele
        gather(key = alls, value = val, -i) %>%
        #   par.locus gives the parent from whom the locus will descend
        mutate(par.locus = gsub('\\_[ab]', '', alls)) %>%
        #   for each locus on each chromosome, pick exactly one parental allele
        group_by(i, par.locus) %>%
        sample_n(size = 1) %>%
        # Remove the unnecessary "parent" column
        select(-alls) %>%
        ungroup() %>%
        mutate(par.locus = gsub('^mom', 'a', par.locus),
               par.locus = gsub('^dad', 'b', par.locus)) %>%
        # Turn this data frame back into "wide" format
        spread(key = par.locus, value = val) %>%
        ungroup() %>%
        # Now, calculate genotype, etc.
        #   for each offspring
        #   (note: to do this, we need to first remove the 'i' 
        #   column in order to calculate g_i)
        select(-i) %>%
        mutate(g_i = apply(., 1, sum) / sqrt(n.loci),
               i = max(popn$i) + 1:nrow(.),
               fem = sample(c(TRUE, FALSE), size = nrow(.), replace = TRUE),
               z_i = rnorm(nrow(.), mean = g_i, sd = sig.e),
               w_i = w.max * exp(-(z_i - theta)^2 / (2*wfitn^2)),
               r_i = rpois(n = nrow(.), lambda = ifelse(fem, 2 * w_i * exp(-alpha * nrow(.)), 0)),
               gen = max(popn$gen) + 1) %>%
        select(i, g_i, z_i, w_i, r_i, fem, gen, all_of(names.array))
    } else {
      # If non-evolving, initialize the population again
      # with size equal to the number of offspring as prescribed by
      # the previous generation's r_i
      offspring = init.sim(a = a, params = params %>% mutate(n.pop0 = sum(popn$r_i))) %>%
        mutate(gen = max(popn$gen) + 1)
    }
    
    return(offspring)
    
  } else {
    return(popn %>% sample_n(size = 0))
  }
  
  
}

##### Testing:
# (again: can not guarantee quality of results)
# (need to run test code above , at least, for this to work)
# # # Try it out.
# propagate.sim(a = c(-1/2, 1/2), params = pars, popn = popn0, evolve = FALSE)
# 
# # Each of the following tests hit the 'if' loop
# # These tests should return an empty data frame.
# 
# # Try it out with an empty data frame.
# propagate.sim(a = c(-1/2, 1/2), params = pars,
#               popn = popn0 %>% sample_n(size = 0))
# # Good.
# 
# # Try it out with only one individual. Should fail.
# propagate.sim(a = c(-1/2, 1/2), params = pars,
#               popn = popn0 %>% sample_n(size = 1))
# # Good.
# 
# # Try it out with two males. Should fail.
# propagate.sim(a = c(-1/2, 1/2), params = pars,
#               popn = popn0[c(1,5),])
# # Good.
# 
# # Try it out with two females. Should fail.
# propagate.sim(a = c(-1/2, 1/2), params = pars,
#               popn = popn0[2:3,])
# # Good.

# Wrapper function for running simulation (using init and propagate functions
# above)
sim = function(a = c(1/2, -1/2), params, init.popn = NULL, evolve = TRUE) {
  
  # N.b., defaults are (1) alleles of equal and opposing effect, (2) no initial
  # population fed in (so init.sim() is called) and (3) populations evolve
  
  # Helpful global parameters:
  
  # How long the simulation should run for
  end.time = params$end.time
  # How many rows the data frame should be initialized with.
  init.row = params$init.row
  # How many loci there are for the allele.
  n.loci = params$n.loci
  # Threshold above which to stop simulation (used in long-sims)
  size.thresh = ifelse('size.thresh' %in% names(params), params$size.thresh, Inf)
  
  
  # A character (string) array for handy indexing
  names.array = paste0(c('a', 'b'), rep(1:n.loci, each = 2))
  
  all.pop = data.frame(i = rep(NA, init.row),
                       g_i = rep(NA, init.row),
                       z_i = rep(NA, init.row),
                       w_i = rep(NA, init.row),
                       r_i = rep(NA, init.row),
                       fem = rep(NA, init.row),
                       gen = rep(NA, init.row)) %>%
    cbind(matrix(NA, 
                 nrow = init.row,
                 ncol = 2 * n.loci) %>%
            data.frame() %>%
            set.names(names.array))
  
  # Initialize the population.
  if (!is.null(init.popn)) { 
    # If the population is passed into the function (i.e., brought from an
    # external environment), its trait value and fitness still need to be
    # determined by information from this environment (the global variables).
    pop0 = init.popn %>%
      mutate(z_i = rnorm(nrow(.), mean = g_i, sd = params$sig.e),
             w_i = params$w.max * exp(-(z_i - params$theta)^2 / (2*params$wfitn^2)),
             r_i = rpois(n = nrow(.), lambda = ifelse(fem, 2 * w_i * exp(-params$alpha * nrow(.)), 0)),
             i = 1:nrow(.),
             gen = 1) %>%
      select(i, g_i, z_i, w_i, r_i, fem, gen, all_of(names.array))
  } else {                   
    pop0 = init.sim(a, params) 
  }
  
  all.pop = dim.add(df = all.pop, 
                    row = init.row,
                    addition = pop0)
  
  prev.gen = pop0
  
  for (time.step in 2:end.time) {
    if(nrow(prev.gen)) {
      if (sum(prev.gen$r_i) < size.thresh) {
        pop = propagate.sim(a = a,
                            params = params,
                            popn = prev.gen,
                            evolve = evolve)
        all.pop = dim.add(df = all.pop,
                          rows = init.row,
                          addition = pop)
        prev.gen = pop
      }
      else return(all.pop %>% filter(!is.na(i)))
    }
  }
  
  # Return population with empty/unused rows filtered out
  return(all.pop %>% filter(!is.na(i)))
  
}

# # More test code below:
# 
# set.seed(12121513)
# 
# sim.test = sim(
#     a = c(-1/2, 1/2),
#     params = data.frame(end.time = 15,
#                         init.row = 1e4,
#                         n.loci = 20,
#                         n.pop0 = 40,
#                         w.max = 1.4,
#                         theta = 1,
#                         wfitn = sqrt(1 / 0.14),
#                         sig.e = 0.5)
# )
# 
# Test to see if this works with inserted initial population.
# set.seed(12121513)
# 
# sim.test = sim(
#     a = c(-1/2, 1/2),
#     params = data.frame(end.time = 15,
#                         init.row = 1e4,
#                         n.loci = 20,
#                         w.max = 1.4,
#                         theta = 2,
#                         wfitn = sqrt(1 / 0.14),
#                         sig.e = 0.5),
#     init.popn = popn0
# )
# 
# sim.test = sim(
#   a = c(-1/2, 1/2),
#   params = data.frame(end.time = 15,
#                       init.row = 1e4,
#                       n.loci = 20,
#                       w.max = 1.4,
#                       theta = 1,
#                       wfitn = sqrt(1 / 0.14),
#                       sig.e = 0.5),
#   init.popn = popn0[1:20,]
# )
# 
# # A test of non-zero alpha
# set.seed(12121513)
# 
# sim.test = sim(
#     a = c(-1/2, 1/2),
#     params = data.frame(end.time = 15,
#                         init.row = 1e4,
#                         n.loci = 20,
#                         n.pop0 = 40,
#                         w.max = 1.4,
#                         theta = 1,
#                         wfitn = sqrt(1 / 0.14),
#                         sig.e = 0.5,
#                         alpha = 0.0035)
# )
# Seems to work!
# With custom insert population and ndd
# sim.test = sim(
#   a = c(-1/2, 1/2),
#   params = data.frame(end.time = 15,
#                       init.row = 1e4,
#                       n.loci = 20,
#                       n.pop0 = 40,
#                       w.max = 1.4,
#                       theta = 1,
#                       wfitn = sqrt(1 / 0.14),
#                       sig.e = 0.5,
#                       alpha = 0.0035),
#   init.popn = popn0[1:20,]
# )

# # Examine output
# sim.test
# 
# # Look at population size over time.
# sim.test %>%
#     group_by(gen) %>%
#     summarise(n = n())
# 
# # Look at allele frequencies over time
# sim.test %>%
#     gather(key = loc.copy, value = val, -c(i, g_i, z_i, w_i, r_i, fem, gen)) %>%
#     mutate(loc = gsub('[ab]', '', loc.copy) %>% as.numeric()) %>%
#     group_by(gen, loc) %>%
#     summarise(pos.freq = mean(val > 0)) %>%
#     #spread(key = loc, val = pos.freq) %>%
#     spread(key = gen, val = pos.freq)
# 
# # Look at genetic variance over time
# sim.test %>%
#   gather(key = loc.copy, value = val, -c(i, g_i, z_i, w_i, r_i, fem, gen)) %>%
#   mutate(loc = gsub('[ab]', '', loc.copy) %>% as.numeric()) %>%
#   group_by(gen, loc) %>%
#   summarise(pos.freq = mean(val > 0),
#             gen.varn = 2 * pos.freq * (1 - pos.freq) / 25) %>%
#   group_by(gen) %>%
#   summarise(gen.var = sum(gen.varn))
# 
