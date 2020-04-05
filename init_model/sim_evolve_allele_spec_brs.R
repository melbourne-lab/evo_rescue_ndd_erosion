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

init.sim = function(a, params) {
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
    w.sqd = params$w.sqd
    # variance of the selection pressure
    
    
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
    #   - b_i:  breeding value; sum of alleles, scaled by sqrt(n)
    #   (scaling by sqrt(n) ensures sums don't explode to infinity)
    #   - i:  a unique identifier for each individual
    #   - fem: assign sex to each individual
    #   - r_i: number of offspring for each individual
    #   Poisson distributed, with mean
    #   W_i = W_max * exp(-k * (b_i - theta)^2)
    init.popn = matrix(sample(a, 2 * n.loci * n.pop0, replace = TRUE), 
           nrow = n.pop0, 
           ncol = 2 * n.loci) %>%
        data.frame() %>%
        set.names(name.array = all_of(names.array)) %>%
        mutate(b_i = apply(., 1, sum) / sqrt(n.loci),
               i = 1:n.pop0,
               fem = sample(c(TRUE, FALSE), size = n.pop0, replace = TRUE),
               w_i = w.max * exp(-(b_i - theta)^2 / w.sqd),
               r_i = rpois(n = n.pop0, lambda = ifelse(fem, 2 * w_i, 0)),
               gen = 1) %>%
        select(i, b_i, w_i, r_i, fem, gen, all_of(names.array))
    
    return(init.popn)
    
}

##### Test of the above:
# 
# pars = data.frame(n.loci = 20, n.pop0 = 40,
#                   w.max = 1.2, theta = 1, w.sqd = 1 / 0.14)
# 
# popn0 = init.sim(a = c(-1/2, 1/2), params = pars)
# popn0
# 
# # Looks good.

propagate.sim = function(a, params, popn) {
    
    n.loci = params$n.loci
    # number of loci determining the genotype
    n.pop0 = params$n.pop0
    # initial population size
    w.max = params$w.max
    # max population size
    theta = params$theta
    # optimal phenotypic value
    w.sqd = params$w.sqd
    # variance of the selection pressure
    
    # A character (string) array for handy indexing
    names.array = paste0(c('a', 'b'), rep(1:n.loci, each = 2))
    
    # NOTE: need a flag to handle an empty array
    # being passed in.
    # (is any(popn$fem) made redundant by the sum?)
    if (any(popn$fem) & any(!popn$fem) & sum(popn$r_i[popn$fem])) {
        
        offspring = cbind(
            # Maternal data frame:
            # Take the female rows in the data frame
            # Remove unnecessary columns (don't need to be inherited)
            # Rename columns to indicate alleles inhereted from mom
            #   NOTE: r_i also included here because we will need it later
            popn %>% 
                filter(fem) %>% 
                select(-c(i, b_i, w_i, fem, gen)) %>%
                set.names(paste(names(.), 
                                ifelse(grepl('^[ab]\\d', names(.)), 'mom', ''), 
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
                set.names(paste(names(.), 'dad', sep = '_'))
            ) %>%
            # For each mating pair, duplicate by the number of offspring
            #   as determined by r_i
            uncount(weight = r_i_) %>%
            # Add a new column for unique ID of each individual
            #   (note - we'll have to remove this later for silly reasons)
            mutate(i = max(popn$i) + 1:nrow(.)) %>%
            # Use some cleverness to segregate alleles:
            #   create a row for each allele
            #   chrom.locus gives the location of the locus
            gather(key = alls, value = val, -i) %>%
            separate(col = alls, into = c('chrom.locus', 'parent'), sep = '_') %>%
            #   for each locus on each chromosome, pick exactly one parental allele
            group_by(i, chrom.locus) %>%
            sample_n(size = 1) %>%
            # Remove the unnecessary "parent" column
            select(-parent) %>%
            # Turn this data frame back into "wide" format
            spread(key = chrom.locus, value = val) %>%
            ungroup() %>%
            # Now, calculate breeding value (genotype?), etc.
            #   for each offspring
            #   (note: to do this, we need to first remove the 'i' 
            #   column in order to calculate b_i)
            select(-i) %>%
            mutate(b_i = apply(., 1, sum) / sqrt(n.loci),
                   i = max(popn$i) + 1:nrow(.),
                   fem = sample(c(TRUE, FALSE), size = nrow(.), replace = TRUE),
                   w_i = w.max * exp(-(b_i - theta)^2 / w.sqd),
                   r_i = rpois(n = nrow(.), lambda = ifelse(fem, 2 * w_i, 0)),
                   gen = max(popn$gen) + 1) %>%
            select(i, b_i, w_i, r_i, fem, gen, all_of(names.array))
            
            return(offspring)
        
    } else {
            return(popn %>% sample_n(size = 0))
    }
        

}

##### Testing:
# # Try it out.
# propagate.sim(a = c(-1/2, 1/2), params = pars, popn = popn0)
# 
# # Try it out with an empty data frame. Should fail.
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

sim = function(a, params) {
    
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
                         b_i = rep(NA, init.row),
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
    pop0 = init.sim(a, params)
    
    all.pop = dim.add(df = all.pop, 
                      row = init.row,
                      addition = pop0)
    
    prev.gen = pop0
    
    for (time.step in 2:end.time) {
        if(nrow(prev.gen)) {
        pop = propagate.sim(a = a,
                            params = params,
                            popn = prev.gen)
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

# More test code below:

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
#                         w.sqd = 1 / 0.14)
# )
# 
# sim.test
# 
# sim.test %>%
#     group_by(gen) %>%
#     summarise(n = n())
# 
# sim.test %>%
#     gather(key = loc.copy, value = val, -c(i, b_i, r_i, fem, gen)) %>%
#     mutate(loc = gsub('[ab]', '', loc.copy) %>% as.numeric()) %>%
#     group_by(gen, loc) %>%
#     summarise(pos.freq = mean(val > 0)) %>%
#     #spread(key = loc, val = pos.freq) %>%
#     spread(key = gen, val = pos.freq)
# 
# sim.test %>%
#     gather(key = loc.copy, value = val, -c(i, b_i, r_i, fem, gen)) %>%
#     mutate(loc = gsub('[ab]', '', loc.copy) %>% as.numeric()) %>%
#     group_by(gen, loc) %>%
#     summarise(g.var = var(val)) %>%
#     #spread(key = loc, val = pos.freq) %>%
#     spread(key = gen, val = g.var)
# 
# sim.test %>%
#     gather(key = loc.copy, value = val, -c(i, b_i, r_i, fem, gen)) %>%
#     mutate(loc = gsub('[ab]', '', loc.copy) %>% as.numeric()) %>%
#     group_by(gen, loc) %>%
#     summarise(n = length(unique(i)),
#               m = mean(val),
#               v = var(val)) %>%
#     print(n = 150)
# 
# sim.test %>%
#     gather(key = loc.copy, value = val, -c(i, b_i, r_i, fem, gen)) %>%
#     mutate(loc = gsub('[ab]', '', loc.copy) %>% as.numeric()) %>%
#     group_by(gen, loc) %>%
#     summarise(n = length(unique(i)),
#               m = mean(val),
#               v = var(val)) %>%
#     group_by(gen) %>%
#     summarise(m = mean(v),
#               n = mean(n)) %>%
#     ggplot() +
#     geom_line(aes(x = gen, y = m))
# 
# 
# object.size(sim.test)
# 
# set.seed(142213)
# 
# sim.test = sim(
#     a = c(-1/2, 1/2),
#     params = data.frame(end.time = 15,
#                         init.row = 1e4,
#                         n.loci = 40, 
#                         n.pop0 = 100,
#                         w.max = 1.4, 
#                         theta = 1, 
#                         w.sqd = 1 / 0.14)
# )
# 
# sim.test %>%
#     gather(key = loc.copy, value = val, -c(i, b_i, w_i, r_i, fem, gen)) %>%
#     mutate(loc = gsub('[ab]', '', loc.copy) %>% as.numeric()) %>%
#     group_by(gen, loc) %>%
#     summarise(n = length(unique(i)),
#               m = mean(val),
#               v = var(val),
#               w = mean(w_i),
#               fixp = all(val > 0),
#               fixn = all(val < 0)) %>%
#     group_by(gen) %>%
#     summarise(m = mean(v),
#               n = mean(n),
#               n.fixp = sum(fixp),
#               n.fixn = sum(fixn)) %>%
#     ggplot() +
#     geom_line(aes(x = gen, y = m))
# 
# sim1 %>%
#     group_by(gen) %>%
#     summarise(w.bar = mean(w_i),
#               w.var = var(w_i),
#               n = n()) %>%
#     ggplot() +
#     geom_line(aes(x = gen, y = w.bar)) +
#     geom_ribbon(aes(x = gen,
#                    ymin = w.bar - 2 * sqrt(w.var),
#                    ymax = w.bar + 2 * sqrt(w.var)),
#                 alpha = 0.2)
