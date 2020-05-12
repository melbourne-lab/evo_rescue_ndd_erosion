# Analyze first time steps in model
# Ran 2000 first-steps for a param combo on Brett's server.

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

# Parameters used

params = data.frame(end.time = 2,
                    init.row = 1e4,
                    n.loci = 25, 
                    n.pop0 = 20,
                    w.max = 2, 
                    theta = 2.6, 
                    wfitn = sqrt(1 / 0.14 / 2),
                    sig.e = 0)

# The simulation data
firsts = read.csv('analyze_model/genotypic changes/first_timers.csv')

# Look at change in allele frequencies 
first.long = firsts %>%
  select(-c(z_i, fem)) %>%
  gather(locus, gval, -c(i, g_i, w_i, r_i, trial, gen)) %>%
  mutate(lcopy = gsub('\\d', '', locus),
         locus = gsub('[ab]', '', locus))

# Allele frequencies - wide form, all gens
p.all.wide = first.long %>%
  group_by(locus, trial, gen) %>%
  summarise(mean.val = mean(gval > 0)) %>%
  mutate(gen = paste0('gen', gen)) %>%
  spread(gen, mean.val)
  
p.all.wide %>%
  ggplot() +
  geom_line(aes(x = gen1, y = gen2, group = trial), size = 0.2)

p.all.wide %>%
  ggplot() +
  geom_point(aes(x = gen1, y = gen2), size = 0.2)

p.all.wide %>% lm(formula = gen2 ~ gen1) %>% summary()
# 0.015629  pm 2 *  0.003291 is the change over time