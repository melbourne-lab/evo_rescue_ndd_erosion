# rm(list = ls())

source('sim_functions.R')

init.sim.bh = function(params, a = c(-1/2, 1/2)) {
  
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
  beta  = ifelse('beta' %in% names(params), params$beta, 0)
  # strength of density dependence (in Beverton-Holt equation)
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
    mutate(
      g_i = apply(., 1, sum) / sqrt(n.loci),
      i = 1:n.pop0,
      fem = sample(c(TRUE, FALSE), size = n.pop0, replace = TRUE),
      z_i = rnorm(n.pop0, mean = g_i, sd = sig.e),
      w_i = w.max * exp(-(z_i - theta)^2 / (2 * wfitn^2)),
      r_i = rpois(n = n.pop0, lambda = ifelse(fem, 2 * w_i / (1 + beta * nrow(.)), 0)),
      gen = 1
    ) %>%
    select(i, g_i, z_i, w_i, r_i, fem, gen, all_of(names.array))
  
  return(init.popn)
  
}

propagate.sim.bh = function(params, a = c(1/2, -1/2), popn, evolve = TRUE) {
  
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
  beta = ifelse('beta' %in% names(params), params$beta, 0)
  # strength of density dependence (in Beverton-Holt formulation)
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
          set.names(
            paste(ifelse(grepl('^[ab]\\d', names(.)), 'mom', ''),
                  names(.), 
                  sep = '_')
          ),
        # Paternal data frame
        # Sample these to get mating pairs, i.e.,
        #   draw from the pool of males once for each female
        # Select the paternal alleles
        # Rename columns to indicate alleles inhereited from dad
        # NOTE: this assumes that each mom mates with only one dad
        popn %>% 
          sample_n(
            size = sum(fem), 
            weight = as.numeric(!fem) / sum(as.numeric(!fem)),
            replace = TRUE
          ) %>%
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
        #   create a row for each locus
        # # gather(key = locs, value = val, -i) %>%
        pivot_longer(-i, names_to = "locs", values_to = "val") %>%
        #   par.locus gives the parent from whom the locus will descend
        mutate(par.locus = gsub('\\_[ab]', '', locs)) %>%
        #   for each locus on each chromosome, pick exactly one parental allele
        group_by(i, par.locus) %>%
        sample_n(size = 1) %>%
        ungroup() %>%
        # Remove the unnecessary allele column
        select(-locs) %>%
        mutate(
          par.locus = gsub('^mom', 'a', par.locus),
          par.locus = gsub('^dad', 'b', par.locus)
        ) %>%
        # Turn this data frame back into "wide" format
        # # spread(key = par.locus, value = val) %>%
        pivot_wider(names_from = par.locus, values_from = val) %>%
        ungroup() %>%
        # Now, calculate genotype, etc.
        #   for each offspring
        #   (note: to do this, we need to first remove the 'i' 
        #   column in order to calculate g_i)
        select(-i) %>%
        mutate(
          g_i = apply(., 1, sum) / sqrt(n.loci),
          i = max(popn$i) + 1:nrow(.),
          fem = sample(c(TRUE, FALSE), size = nrow(.), replace = TRUE),
          z_i = rnorm(nrow(.), mean = g_i, sd = sig.e),
          w_i = w.max * exp(-(z_i - theta)^2 / (2*wfitn^2)),
          r_i = rpois(n = nrow(.), lambda = ifelse(fem, 2 * w_i / (1 + beta * nrow(.)), 0)),
          gen = max(popn$gen) + 1
        ) %>%
        select(i, g_i, z_i, w_i, r_i, fem, gen, all_of(names.array))
    } else {
      # If non-evolving, initialize the population again
      # with size equal to the number of offspring as prescribed by
      # the previous generation's r_i
      offspring = init.sim.bh(a = a, params = params %>% mutate(n.pop0 = sum(popn$r_i))) %>%
        mutate(gen = max(popn$gen) + 1)
    }
    
    return(offspring)
    
  } else {
    return(popn %>% sample_n(size = 0))
  }
  
  
}


sim.bh = function(params, a = c(1/2, -1/2), init.popn = NULL, evolve = TRUE) {
  
  # N.b., defaults are (1) alleles of equal and opposing effect, (2) no initial
  # population fed in (so init.sim() is called) and (3) populations evolve
  
  # Helpful global parameters:
  
  # How long the simulation should run for
  end.time = params$end.time
  # How many rows the data frame should be initialized with.
  init.row = params$init.row
  # How many loci there are
  n.loci   = params$n.loci
  # Threshold above which to stop simulation (used in long-sims)
  size.thresh = ifelse('size.thresh' %in% names(params), params$size.thresh, Inf)
  
  
  # A character (string) array for handy indexing
  names.array = paste0(c('a', 'b'), rep(1:n.loci, each = 2))
  
  all.pop = data.frame(
    i = rep(NA, init.row),
    g_i = rep(NA, init.row),
    z_i = rep(NA, init.row),
    w_i = rep(NA, init.row),
    r_i = rep(NA, init.row),
    fem = rep(NA, init.row),
    gen = rep(NA, init.row)
  ) %>%
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
    pop0 = data.frame(
      g_i = init.popn %>%
        select(all_of(names.array)) %>%
        (function(x) apply(x, 1, sum) / sqrt(n.loci))
    ) %>%
      mutate(
        fem = init.popn$fem,
        z_i = rnorm(nrow(.), mean = g_i, sd = params$sig.e),
        w_i = params$w.max * exp(-(z_i - params$theta)^2 / (2*params$wfitn^2)),
        r_i = rpois(n = nrow(.), lambda = ifelse(fem, 2 * w_i * (1 + params$beta * nrow(.)), 0)),
        i = 1:nrow(.),
        gen = 1
      ) %>%
      cbind(init.popn %>% select(all_of(names.array))) %>%
      select(i, g_i, z_i, w_i, r_i, fem, gen, all_of(names.array))
  } else {                   
    pop0 = init.sim.bh(params, a) 
  }
  
  all.pop = dim.add(df = all.pop, 
                    row = init.row,
                    addition = pop0)
  
  prev.gen = pop0
  
  for (time.step in 2:end.time) {
    if(nrow(prev.gen)) {
      if (sum(prev.gen$r_i) < size.thresh) {
        pop = propagate.sim.bh(params = params,
                               a = a,
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
