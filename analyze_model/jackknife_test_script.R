# just kinda want to see what the effect of one individual on the whole population is
# using a crude jackknifing approach here
# (originally run on the server but git auth issues make that a headache)
# need to think about this more
# SN - 5-12 April 2022

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

replace.i = function(i.rm, popn, sex.same = FALSE) {
  new.i = max(popn$i) + 1
  i.fem = popn$fem[popn$i %in% i.rm]
  popn  = popn %>% filter(!i %in% i.rm)
  
  if (!sex.same) { repl = popn %>% sample_n(1)
  } else {         repl = popn %>% filter(fem %in% i.fem) %>% sample_n(1) }
  
  popn.out = rbind(popn, repl %>% mutate(i = new.i))
  
  return(popn.out)
}

source('base_model_source/sim_functions.R')

pop.size = 100
n.trials = 100

pars = data.frame(
  end.time = 30,
  n.loci = 25,
  n.pop0 = pop.size,
  w.max  = 2, 
  theta  = 2.8,
  wfitn  = sqrt(1 / 0.14 / 2),
  sig.e  = sqrt(0.5),
  pos.p  = 0.5,
  alpha  = 0.0035,
  init.row = 1e5
)

set.seed(209)

source.pop = init.sim(params = pars)

full.list = vector('list', n.trials)

set.seed(490)

for (j in 1:n.trials) {
  
  sim.output = sim(params = pars, init.popn = source.pop)
  
  demo.summ = sim.output %>%
    group_by(gen) %>%
    summarise(n = n(),
              gbar = mean(g_i),
              zbar = mean(z_i),
              wbar = mean(w_i))
  
  gene.summ = sim.output %>%
    select(-c(g_i, z_i, w_i, r_i, fem)) %>%
    gather(key = loc.copy, value = val, -c(i, gen)) %>%
    mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
    group_by(i, gen, locus) %>%
    mutate(heter = !sum(val)) %>%
    group_by(gen, locus) %>%
    summarise(p = mean(val > 0)) %>%
    group_by(gen) %>%
    summarise(p.fix.pos = mean(p == 1),
              p.fix.neg = mean(p == 0),
              v = sum(2 * p * (1 - p)) / pars$n.loci[1])
  
  full.list[[j]] = cbind(demo.summ, gene.summ %>% select(-gen)) %>%
    mutate(trial.no = j)
  
  print(j)
  
}

full.data = do.call(what = rbind, full.list)

full.data %>%
  select(gen, trial.no, n) %>%
  rbind(expand.grid(gen = 1:pars$end.time,
                    trial.no = 1:n.trials,
                    n = 0)) %>%
  group_by(gen, trial.no) %>%
  summarise(n = sum(n)) %>%
  group_by(gen) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    nn = n()) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = nbar
    )
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / nn),
      ymax = nbar + 2 * sqrt(nvar / nn)
    ),
    alpha = 0.1
  ) +
  geom_line(
    data = full.data,
    aes(
      y = n,
      group = trial.no
    ),
    size = 0.1
  ) +
  scale_y_log10()

### Jackknifing

jack.trials = expand.grid(i = 1:pars$n.pop0, trial = 1:n.trials)

jack.list = vector('list', nrow(jack.trials))

set.seed(10001)

for (j in 1:nrow(jack.trials)) {
  
  # i.removed = sim.trials$i[j]
  # trial.no  = sim.trials$trial[j]
  
  sim.output = sim(params = pars,
                   init.popn = replace.i(i.rm = jack.trials$i[j],
                                         popn = source.pop,
                                         sex.same = TRUE))
  
  demo.summ = sim.output %>%
    group_by(gen) %>%
    summarise(n = n(),
              gbar = mean(g_i),
              zbar = mean(z_i),
              wbar = mean(w_i))
  
  gene.summ = sim.output %>%
    select(-c(g_i, z_i, w_i, r_i, fem)) %>%
    gather(key = loc.copy, value = val, -c(i, gen)) %>%
    mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
    group_by(i, gen, locus) %>%
    mutate(heter = !sum(val)) %>%
    group_by(gen, locus) %>%
    summarise(p = mean(val > 0)) %>%
    group_by(gen) %>%
    summarise(p.fix.pos = mean(p == 1),
              p.fix.neg = mean(p == 0),
              v = sum(2 * p * (1 - p)) / pars$n.loci[1])
  
  jack.list[[j]] = cbind(demo.summ, gene.summ %>% select(-gen)) %>%
    mutate(i.removed = jack.trials$i[j],
           trial.no  = jack.trials$trial[j])
  
  print(paste(j, 'of', nrow(jack.trials)))
    
}

# `summarise()` has grouped output by 'gen'. You can override using the `.groups` argument.
# [1] "9864 of 10000"
# Error in select(., -alls) : bad generic call environment
### thonkas...

jack.all = merge(
  x = do.call(what = rbind, jack.list),
  y = jack.trials,
  by.x = c('i.removed', 'trial.no'),
  by.y = c('i', 'trial')
) 

jack.all %>%
  group_by(i.removed, gen) %>%
  summarise(nbar = mean(n)) %>%
  ggplot(aes(x = gen, y = nbar, group = i.removed)) +
  geom_line()

# need to add zeros

jack.n = jack.all %>%
  select(i.removed, trial.no, n, gen) %>%
  rbind(expand.grid(i.removed = 1:pars$n.pop0,
                    trial.no  = 1:n.trials,
                    n = 0,
                    gen = 1:pars$end.time)) %>%
  group_by(i.removed, trial.no, gen) %>%
  summarise(n = sum(n)) %>%
  filter(any(n > 0))

jack.n %>%
  group_by(i.removed, gen) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            nn = n()) %>%
  merge(y = source.pop %>% select(i, g_i),
        by.x = 'i.removed', by.y = 'i') %>%
  ggplot(aes(x = gen, y = nbar)) +
  geom_line(aes(colour = g_i, group = i.removed)) +
  scale_colour_viridis_c()
  
jack.n %>%
  group_by(i.removed, trial.no) %>%
  summarise(extinct = any(!n)) %>%
  group_by(i.removed) %>%
  summarise(pext = mean(extinct)) %>%
  merge(y = source.pop %>% select(i, g_i),
        by.x = 'i.removed', by.y = 'i') %>%
  ggplot(aes(x = g_i, y = pext)) +
  geom_point()

