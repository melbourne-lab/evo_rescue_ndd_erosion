library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

# Auxiliary functions

set.names = function(df, name.array) {
  names(df) = name.array
  df
}

dim.add = function(df, rows, addition) {
  if (sum(is.na(df$i)) < nrow(addition)) {
    df = df %>%
      rbind(df %>%
              sample_n(size = rows) %>%
              mutate_all(function(x) NA))
  }
  if (nrow(addition)) {
    df[which.max(is.na(df$i)) + 0:(nrow(addition)-1),] = addition
  }
  return(df)
}

init.sim.overlap = function(a = c(1/2, -1/2), params) {
  # Inputs:
  # a - an array of length two (bi-allelic model)
  # each element is a contribution to the genotype
  # in simplest case, this is -1/2, 1/2
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
           r_i = rpois(n = n.pop0, lambda = ifelse(fem, w_i * exp(-alpha * n.pop0), 0)),
           age = sample(1:2, size = n.pop0, replace = TRUE),
           gen = 1) %>%
    select(i, g_i, z_i, w_i, r_i, fem, age, gen, all_of(names.array))
  
  return(init.popn)
  
}

# #### Test of the above:
# 
# pars = data.frame(n.loci = 20, n.pop0 = 40,
#                   w.max = 1.2, theta = 2.6,
#                   wfitn = sqrt(1 / 0.14),
#                   sig.e = 0.5)
# 
# popn0 = init.sim.overlap(a = c(-1/2, 1/2), params = pars)
# popn0

# Looks good.

# # Trying a sim with a custom allele frequency:
# 
# pars = data.frame(n.loci = 20, n.pop0 = 40,
#                   w.max = 1.2, theta = 2.6,
#                   wfitn = sqrt(1 / 0.14),
#                   sig.e = 0.5, pos.p = 0.5)
# 
# popn0 = init.sim(params = pars)
# popn0

propagate.sim.overlap = function(a = c(1/2, -1/2), params, popn, evolve = TRUE) {
  
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
  #   - any individuals age 1 (to survive to next generation)
  if (with(popn, any(fem) & any(!fem) & sum(r_i))) {
    
    if (evolve) {
      
      offspring = cbind(
        # Maternal data frame:
        # Take the female rows
        # Remove unnecessary columns (don't need to be inherited)
        # Rename columns to indicate alleles inhereted from mom
        #   NOTE: r_i also included here because we will need it later
        popn %>% 
          filter(fem) %>% 
          select(-c(i, g_i, w_i, z_i, fem, age, gen)) %>%
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
        # Now, calculate breeding value (genotype?), etc.
        #   for each offspring
        #   (note: to do this, we need to first remove the 'i' 
        #   column in order to calculate g_i)
        select(-i) %>%
        mutate(g_i = apply(., 1, sum) / sqrt(n.loci),
               i = max(popn$i) + 1:nrow(.),
               fem = sample(c(TRUE, FALSE), size = nrow(.), replace = TRUE),
               z_i = rnorm(nrow(.), mean = g_i, sd = sig.e),
               w_i = w.max * exp(-(z_i - theta)^2 / (2*wfitn^2)),
               r_i = rpois(n = nrow(.), lambda = ifelse(fem, w_i * exp(-alpha * nrow(.)), 0)),
               age = 1,
               gen = max(popn$gen) + 1) %>%
        select(i, g_i, z_i, w_i, r_i, fem, age, gen, all_of(names.array)) 
      
      if (any(popn$age < 2)) {
        
        next.gen = offspring %>%
          # Add one year old parents from above
          #   take prev. generaton, filter out 
          #   update age from 1 to 2 and update generation
          #   recalculate number of offspring
          rbind(
            popn %>%
              filter(age < 2) %>%
              mutate(age = age + 1,
                     gen = gen + 1) %>%
              mutate(r_i = rpois(n = nrow(.), lambda = ifelse(fem, w_i * exp(-alpha * nrow(.)), 0)))
          )
        
        return(next.gen)
      } else { return(offspring) }
      
    } else {
      # If non-evolving, initialize the population again
      # with size equal to the number of offspring as prescribed by
      # the previous generation's r_i
      offspring = init.sim.overlap(a = a, params = params %>% mutate(n.pop0 = sum(popn$r_i))) %>%
        mutate(gen = max(popn$gen) + 1)
    }
    
    return(offspring)
    
  } 
  
  else {
    return(popn %>% sample_n(size = 0))
  }
  
}

# ##### Testing:
# # # # Try it out.
# popn1 = propagate.sim.overlap(a = c(-1/2, 1/2), params = pars, popn = popn0, evolve = TRUE)
# popn2 = propagate.sim.overlap(a = c(-1/2, 1/2), params = pars, popn = popn1, evolve = TRUE)
# 
# # Each of the following tests hit the 'if' loop
# # These tests should return an empty data frame.
# 
# # Try it out with an empty data frame.
# propagate.sim.overlap(a = c(-1/2, 1/2), params = pars,
#                       popn = popn0 %>% sample_n(size = 0))
# # Good.
# 
# # Try it out with only one individual. Should fail.
# propagate.sim.overlap(a = c(-1/2, 1/2), params = pars,
#                       popn = popn0 %>% sample_n(size = 1))
# # Good.
# 
# # Try it out with two males. Should fail.
# propagate.sim.overlap(a = c(-1/2, 1/2), params = pars,
#                       popn = popn0[c(1,5),])
# # Good.
# 
# # Try it out with two females. Should fail.
# propagate.sim.overlap(a = c(-1/2, 1/2), params = pars,
#                       popn = popn0[3:4,])
# # Good.
# 
# # Try feeding in only individuals at age 1 (should still work)
# propagate.sim.overlap(a = c(-1/2, 1/2), params = pars,
#                       popn = popn0 %>% filter(age < 2))
# # Good.
# 
# # Try feeding in only individuals at age 2 (should also work)
# propagate.sim.overlap(a = c(-1/2, 1/2), params = pars,
#                       popn = popn0 %>% filter(age > 1))

sim.overlap = function(a = c(1/2, -1/2), params, init.popn = NULL, evolve = TRUE) {
  
  # Helpful global parameters.
  
  # How long the simulation should run for
  end.time = params$end.time
  # How many rows the data frame should be initialized with.
  init.row = params$init.row
  # How many loci there are for the allele.
  n.loci = params$n.loci
  
  # A character (string) array for handy indexing
  names.array = paste0(c('a', 'b'), rep(1:n.loci, each = 2))
  
  all.pop = data.frame(i = rep(NA, init.row),
                       g_i = rep(NA, init.row),
                       z_i = rep(NA, init.row),
                       w_i = rep(NA, init.row),
                       r_i = rep(NA, init.row),
                       fem = rep(NA, init.row),
                       age = rep(NA, init.row),
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
             r_i = rpois(n = nrow(.), lambda = ifelse(fem, w_i * exp(-params$alpha * nrow(.)), 0)),
             i = 1:nrow(.),
             gen = 1) %>%
      select(i, g_i, z_i, w_i, r_i, fem, age, gen, all_of(names.array))
  } else {                   
    pop0 = init.sim.overlap(a, params) 
  }
  
  all.pop = dim.add(df = all.pop, 
                    row = init.row,
                    addition = pop0)
  
  prev.gen = pop0
  
  for (time.step in 2:end.time) {
    if(nrow(prev.gen)) {
      pop = propagate.sim.overlap(
        a = a,
        params = params,
        popn = prev.gen,
        evolve = evolve
      )
      all.pop = dim.add(df = all.pop,
                        rows = init.row,
                        addition = pop)
      prev.gen = pop
    }
  }
  
  all.pop = all.pop %>% filter(!is.na(i))
  
  return(all.pop)
  
}

unroller = function(sim.list) {
  sim.list %>%
    do.call(what = rbind) %>%
    mutate(trial = cumsum(c(1, as.numeric(diff(gen) < 0))))
}

# # More test code below:
# 
# set.seed(12121513)
# 
# sim.test = sim.propagate(
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
