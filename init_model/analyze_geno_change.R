# A script looking at allele-specific model
# Why and how do variance measures change over time?
# Follow gomulkiewicz and holt?


library(ggplot2)
library(dplyr)
library(tidyr)

k.sim.summary = read.csv('model_building/eg_1k_n0_20.csv')

w2 = (1/0.14) / 2
s2 = 0.5

k = w2 / (w2 + s2)

ks = k.sim.summary %>%
    mutate(g.load = 2.6 - b.bar)

ggplot(ks, aes(x = gen)) +
    geom_line(aes(y = g.load)) +
    geom_line(aes(y = 2.6 * k^(gen-1)), colour = 'purple')

ks$g.load[2] / 2.6 # suggests that k = 0.946
ks$g.load[3] / ks$g.load[2] # 

ks$g.load %>% log() %>% diff() %>% plot(type = 'l', ylim = c(-0.07, 0))

ks %>%
    mutate(lgdiff = g.load %>% log() %>% (function(x) c(diff(x), 0)),
           kestmm = w2 / (w2 + v.bar),
           kestmu = kestmm + (exp(lgdiff)[1] - kestmm[1])) %>%
    filter(gen < 15) %>%
    ggplot(aes(x = gen)) +
    geom_line(aes(y = exp(lgdiff))) +
    geom_line(aes(y = kestmm), colour = 'blue') +
    geom_line(aes(y = kestmu), colour = 'purple')
#    geom_smooth(aes(y = exp(lgdiff)))

# Look at one trial. Look at breeding values versus other.

# If s^2 i is shrinking, then w^2 / (w^2 + s^2) will grow
# This makes sense.

# Set seed for reproducibility
set.seed(40904)

# Run one simulation
test.sim = sim(a = c(-1/2, 1/2),
               params = data.frame(end.time = 15,      # Length of simulatin
                                   init.row = 1e4,     # Number of rows to initialize data
                                   n.loci = 25,        # Number of loci
                                   n.pop0 = 20,        # Initial population size
                                   w.max = 2,          # Max fitness of an individual
                                   theta = 2.6,        # Distance from phenotypic optimum
                                   w.sqd = 1 / 0.14))  # Selection strength

test.sim.alleles = test.sim %>%
    gather(key = chr.locus, val = val, -c(i, b_i, w_i, r_i, fem, gen)) %>%
    mutate(locus = gsub('^[ab]', 'locus', chr.locus)) %>%
    group_by(gen) %>%
    summarise(a.var = var(val))

test.sim.genvars = test.sim %>%
    group_by(gen) %>%
    summarise(g.var = var(b_i),
              n = n()) %>%
    merge(y = test.sim.alleles)

test.sim.genvars %>%
    ggplot(aes(x = gen)) +
    geom_line(aes(y = a.var), colour = 'darkgreen') +
    geom_line(aes(y = g.var), colour = 'maroon') +
    labs(y = 'Variance measure')
# THese are obviously measuring wildly different things at wildly different scales... but what?

# Is genetic value variance really changing that much over time???

test.sim %>% 
    filter(gen %in% c(1,15)) %>% 
    ggplot() + 
    geom_density(aes(x = b_i, group = gen))

# whoah yes

# Okay... time to perhaps actually measure "breeding values" here in a rigorous way.

# Set seed for reproducibility
set.seed(40904)

# Run one simulation
test.sim = sim(a = c(-1/2, 1/2),
               params = data.frame(end.time = 15,      # Length of simulatin
                                   init.row = 1e4,     # Number of rows to initialize data
                                   n.loci = 25,        # Number of loci
                                   n.pop0 = 20,        # Initial population size
                                   w.max = 2,          # Max fitness of an individual
                                   theta = 2.6,        # Distance from phenotypic optimum
                                   w.sqd = 1 / 0.14))  # Selection strength

test.sim

test.sim.vars = test.sim %>%
    select(-c(w_i, r_i, fem)) %>%
    gather(key = allele.copy, value = gen.val, -c(i, b_i, gen)) %>%
    mutate(allele = allele.copy %>% gsub(pattern = '^[ab]', replacement = '') %>% as.numeric()) %>%
    group_by(gen, allele) %>%
    summarise(p.pos = mean(gen.val > 0),
              b.var = var(b_i)) %>%
    mutate(sigma.a = 2 * p.pos * (1 - p.pos) / 25) %>%
    group_by(gen) %>%
    summarise(sigma.a = sum(sigma.a),
              sigma.b = b.var[1]) 

# test.sim.vars %>%
#     ggplot() +
#     geom_line(aes(x = gen, y = sigma.a), colour = 'green4') +
#     geom_line(aes(x = gen, y = sigma.b), colour = 'purple')

test.sim %>%
    group_by(gen) %>%
    summarise(w.bar = mean(w_i),
              z.bar = mean(b_i),
              z.var = var(b_i)) %>%
    merge(y = test.sim.vars) %>%
    mutate(z.bar = z.bar - 2.6) %>%
    mutate(k_t = cumprod((1/0.14/2) / (sigma.a + (1/0.14/2)))) %>%
    mutate(z.bar.pred = z.bar[1] * c(1, k_t[0:14])) %>%
    ggplot(aes(x = gen)) +
    geom_line(aes(y = z.bar)) +
    geom_ribbon(aes(ymin = z.bar - 2 * sqrt(z.var),
                    ymax = z.bar + 2 * sqrt(z.var)),
                alpha = 0.2) +
    geom_line(aes(y = z.bar.pred), colour = 'purple')

# Damn bro wtf? How did genetic value... increase?

test.sim %>%
    group_by(gen) %>%
    summarise(b.var = var(b_i)) %>%
    ggplot() +
    geom_line(aes(x = gen, y = b.var))

test.sim %>%
    ggplot(aes(x = gen, y = b_i)) +
    geom_point(position = position_jitter(width = 0.125)) +
    geom_smooth()

## Now look at all genetic variance changes in 1000K sims
# (summarised on zerver)

all.sigma.a = read.csv('model_building/allelic_all_sigma_a.csv')

all.sigma.a %>%
    filter(trial %in% 1:50) %>%
    group_by(trial) %>%
    mutate(extinct = max(gen) < 15) %>%
    ggplot(aes(x = gen)) +
    geom_line(aes(y = sigma.a, group = trial, colour = extinct),
              size = 0.4) +
    scale_color_manual(values = c('black', 'red'))

all.sigma.a %>%
    filter(trial %in% 1:50) %>%
    group_by(trial) %>%
    mutate(extinct = max(gen) < 15) %>%
    filter(!extinct) %>%
    ggplot(aes(x = gen)) +
    geom_line(aes(y = sigma.a, group = trial))

all.sigma.a %>%
    group_by(trial) %>%
    mutate(extinct = max(gen) < 15) %>%
    group_by(extinct, gen) %>%
    summarise(sig.bar = mean(sigma.a),
              sig.var = var(sigma.a),
              n = n()) %>%
    ggplot() +
    geom_line(aes(x = gen, y = sig.bar, 
                  group = extinct,
                  colour = extinct)) +
    geom_ribbon(aes(x = gen,
                    ymin = sig.bar - 2 * sqrt(sig.var / n),
                    ymax = sig.bar + 2 * sqrt(sig.var / n),
                    fill = extinct,
                    group = extinct),
                alpha = 0.2) +
    scale_color_manual(values = c('black', 'red')) +
    scale_fill_manual(values = c('black', 'red')) 

all.sigma.a %>%
    group_by(gen) %>%
    summarise(sig.bar = mean(sigma.a),
              sig.var = var(sigma.a),
              n = n()) %>%
    ggplot() +
    geom_line(aes(x = gen, y = sig.bar)) +
    geom_point(aes(x = gen, y = sig.bar,
                   size = n / max(n))) +
    geom_ribbon(aes(x = gen,
                    ymin = sig.bar - 2 * sqrt(sig.var / n),
                    ymax = sig.bar + 2 * sqrt(sig.var / n)),
                alpha = 0.2)

nrow(all.sigma.a)

all.sigma.a %>%
    group_by(gen) %>%
    summarise(sig.bar = mean(sigma.a),
              sig.var = var(sigma.a),
              n = n()) %>%
    ungroup() %>%
    mutate(ez_t = 2.6 * c(1, cumprod((1/0.14/2) / (sig.bar[-15] + (1/0.14/2)))),
           ez_t_min = 2.6 * c(1, cumprod((1/0.14/2) / 
                                             (sig.bar[-15] - 2*sqrt(sig.var[-15] / n[-15])
                                              + (1/0.14/2)))),
           ez_t_max = 2.6 * c(1, cumprod((1/0.14/2) /
                                             (sig.bar[-15] + 2*sqrt(sig.var[-15] / n[-15]) 
                                              + (1/0.14/2))))) %>%
    merge(y = k.sim.summary %>% select(-n.trials), by = 'gen') %>%
    ggplot(aes(x = gen)) +
    geom_ribbon(aes(ymin = (2.6 - b.bar) - 2 * sqrt(b.var / n),
                    ymax = (2.6 - b.bar) + 2 * sqrt(b.var / n)),
                alpha = 0.2, fill = 'purple') +
    geom_ribbon(aes(ymin = ez_t_min,
                    ymax = ez_t_max),
                alpha = 0.2) +
    geom_line(aes(y = 2.6 - b.bar), colour = 'purple') +
    geom_point(aes(y = 2.6 - b.bar, colour = n / 1000)) +
    geom_line(aes(y = ez_t)) +
    #scale_color_gradient(low = 'red', high = 'black')
    scale_color_viridis_c(option = 'B') +
    scale_y_log10()
    

k.sim.summary 

# For each female, look for each locus at the number of offspring
# as a function of phenotype.

locus.ws = test.sim %>%
    filter(fem) %>%
    select(-c(r_i, fem)) %>%
    gather(key = 'allele.copy', value = 'val', -c(i, b_i, w_i, gen)) %>%
    mutate(locus = paste0('l', gsub('^[ab]', '', allele.copy))) %>%
    group_by(i, locus, gen) %>%
    summarise(w_i = w_i[1],
              loc.val = sum(val > 0)) %>%
    group_by(gen, locus, loc.val) %>%
    summarise(w_bar = mean(w_i))

locus.ws %>%
    filter(gen %in% 1) %>%
    ggplot() +
    geom_point(aes(x = loc.val, y = w_bar),
               position = position_jitter(width = 0.02),
               alpha = 0.5)
locus.ws %>%
    ggplot() +
    geom_point(aes(x = loc.val, y = w_bar, colour = gen),
               position = position_jitter(width = 0.08),
               alpha = 0.75)

locus.ws %>%
    group_by(gen) %>%
    mutate(cwbar = w_bar - mean(w_bar)) %>%
    ggplot() +
    geom_point(aes(x = loc.val, y = cwbar, colour = gen),
               position = position_jitter(width = 0.2))

locus.ws %>%
    group_by(gen, loc.val) %>%
    summarise(w_bar = mean(w_bar)) %>%
    group_by(gen) %>%
    mutate(cwbar = w_bar - mean(w_bar)) %>%
    ggplot() +
    geom_line(aes(x = loc.val, y = cwbar, group = gen),
              alpha = 0.5) +
    geom_point(aes(x = loc.val, y = cwbar, colour = gen),
               size = 2) +
    scale_color_viridis_c(option = 'B')


locus.ws %>%
    lm(formula = w_bar ~ loc.val + as.factor(gen)) %>%
    summary()

locus.ws %>%
    lm(formula = w_bar ~ loc.val * as.factor(gen)) %>%
    summary()
