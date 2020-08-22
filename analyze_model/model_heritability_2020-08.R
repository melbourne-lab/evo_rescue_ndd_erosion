library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

source('base_model_source/sim_functions_pedigreed.R')

parpar = data.frame(n.pop0 = 100,
                    end.time = 10,
                    init.row = 1e4,
                    n.loci = 25,
                    w.max = 1.1,
                    theta = 0,
                    wfitn = sqrt(1 / 0.14 / 2),
                    sig.e = 0.5,
                    pos.p = 0.5,
                    alpha = 0.000)

herit.test = sim.par.id(a = c(1/2, -1/2),
                        params = parpar)

head(herit.test)

g12 = herit.test %>%
  select(i, g_i, z_i, gen, mom_i, dad_i) %>%
  filter(gen %in% 1:2) #%>%
  merge(x = . %>% filter(gen %in% 1),
        y = . %>% filter(gen %in% 2) %>% select(mom_i, g_i, z_i), 
        by.x = 'i', by.y = 'mom_i')

g12
  
midpars.g1 = merge(x = g12 %>% filter(gen %in% 2),
      y = g12 %>% filter(gen %in% 1) %>% select(i, g_i, z_i),
      # x is offspring, y is parent
      by.x = 'mom_i', by.y = 'i', suffixes = c('', '_mom')) %>%
  merge(y = g12 %>% filter(gen %in% 1) %>% select(i, g_i, z_i),
        by.x = 'dad_i', by.y = 'i', suffixes = c('', '_dad')) %>%
  mutate(z_midp = (z_i_mom + z_i_dad) / 2,
         g_midp = (g_i_mom + g_i_dad) / 2)

midpars.g1 %>%
  ggplot(aes(x = z_midp, y = z_i)) +
  geom_point() +
  stat_smooth(method = 'lm')

midpars.g1 %>%
  ggplot(aes(x = g_midp, y = g_i)) +
  geom_point() +
  stat_smooth(method = 'lm')

g_midp_reg = lm(g_midp ~ g_i, midpars.g1)
summary(g_midp_reg)

gh_t = 2:parpar$end.time
zh_t = 2:parpar$end.time

for (timestep in 1:(parpar$end.time-1)) {
  gens = herit.test %>%
    select(i, g_i, z_i, gen, mom_i, dad_i) %>%
    filter(gen %in% (timestep + 0:1))
  
  midpars = merge(x = gens %>% filter(gen %in% (timestep + 1)),
                  y = gens %>% filter(gen %in% (timestep)) %>% select(i, g_i, z_i),
                  # x is offspring, y is parent
                  by.x = 'mom_i', by.y = 'i', suffixes = c('', '_mom')) %>%
    merge(y = gens %>% filter(gen %in% (timestep)) %>% select(i, g_i, z_i),
          by.x = 'dad_i', by.y = 'i', suffixes = c('', '_dad')) %>%
    mutate(z_midp = (z_i_mom + z_i_dad) / 2,
           g_midp = (g_i_mom + g_i_dad) / 2)
  
  gh_t[timestep] = summary(lm(g_midp ~ g_i, midpars))$coefficients[,'Estimate']['g_i']
  zh_t[timestep] = summary(lm(z_midp ~ z_i, midpars))$coefficients[,'Estimate']['z_i']
}

gh_t

mean(gh_t)
plot(gh_t, type = 'l')

zh_t
mean(zh_t)

herit.test %>%
  select(-c(i, g_i, z_i, w_i, r_i, fem, mom_i, dad_i)) %>%
  gather(key = loc.copy, value = val, -c(gen)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(gen, locus) %>%
  summarise(p = mean(val > 0)) %>%
  group_by(gen) %>%
  summarise(v = sum(2 * p * (1 - p)) / parpar$n.loci[1]) %>%
  plot(type = 'l')

# Heritability declines, as does genetic variation
# I wonder why genetic variation declines?

herit.test %>%
  select(-c(i, g_i, z_i, w_i, r_i, fem, mom_i, dad_i)) %>%
  gather(key = loc.copy, value = val, -c(gen)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(gen, locus) %>%
  summarise(p = mean(val > 0)) %>%
  ggplot() + 
  geom_line(
    aes(x = gen,
        y = p,
        group = locus)
  )

# Oh shit... genetic variaton declines naturally just due to drift
# Fuck fuck fuck fuck fuck... ahhh fuck. Is this bad?
# Do I need mutation?

# Without this genetic variation will decline to zero over long timescales.

# What about straight up ratio of genotypic variance to phenotypic variance.

herit.test %>%
  group_by(gen) %>%
  summarise(h_t = var(g_i) / var(z_i))
# see... these numbers are very different!
# why is that?

herit.test %>%
  group_by(gen) %>%
  summarise(h_t = var(g_i) / var(z_i)) %>%
  plot(type = 'l')
