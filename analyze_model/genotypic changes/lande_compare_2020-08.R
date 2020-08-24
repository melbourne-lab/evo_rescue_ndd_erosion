# Another script where I compare phenotypic predictions with expectations in Lande, 1976.
# August 20 2020

set.seed(4199)

parpar = data.frame(n.pop0 = 100,
                    end.time = 15,
                    init.row = 1e4,
                    n.loci = 25,
                    w.max = 2,
                    theta = 3,
                    wfitn = sqrt(1 / 0.14 / 2),
                    sig.e = 0.5,
                    pos.p = 0.5,
                    alpha = 0.000)

testsim = sim(a = c(1/2, -1/2), params = parpar)

table(testsim$gen)
table(testsim$gen) %>% plot()

tsgen1 = testsim %>% filter(gen %in% 1)

tsgen1 %>% mutate(zwt = z_i * w_i / mean(w_i)) %>% summarise(zbar = mean(z_i), zwbart = parpar$theta - mean(zwt), sig2 = var(z_i))
parpar$wfitn^2 / (parpar$wfitn^2 + 0.908) * (parpar$theta - 0.07168)

# according to this, the sim model gets zwt right

(
tsgen1.summ = testsim %>%
  filter(gen %in% 1:2) %>%
  group_by(gen) %>%
  mutate(d_i = parpar$theta - z_i) %>%
  summarise(dbar = mean(d_i),
            sig2 = var(g_i),
            dwbar = mean(d_i * w_i / mean(w_i)),
            predzbar = dbar * with(parpar, (wfitn^2 + sig.e^2) / (wfitn^2 + sig.e^2 + sig2)))
)
# but, here, the model is off for zbar predicted by Lande 

(
  tsgen.summ = testsim %>%
    group_by(gen) %>%
    mutate(d_i = parpar$theta - z_i) %>%
    summarise(dbar = mean(d_i),
              # sig2 = var(g_i),
              dwbar = mean(d_i * w_i / mean(w_i)),
              preddw = dbar * with(parpar, wfitn^2 / (wfitn^2 + var(z_i))),
              preddbar = dbar * with(parpar, (wfitn^2 + sig.e^2) / (wfitn^2 + sig.e^2 + var(g_i)))) #%>%
    # mutate(prS = (dbar - dwbar),
    #        prh = c(-diff(dbar), 0) / prS) %>%
)

tsgen.summ$prh %>% mean()

plot(tsgen.summ$preddbar[1:14], tsgen.summ$dbar[2:15])
abline(0,1)
# as we've seen, prediction is biased
plot(dwbar ~ preddw, tsgen.summ)
abline(0,1)
# but the dwbar is correct... this means Lande's z_w(t) is correct

# Maybe heritability in our model is not actual h^2?
# Test heritability in model tomorrow.
# Also wonder if the sex does anything - should increase variance but not bias results...

### let's try this with pedigreed sims

rm(list = ls())

source("base_model_source/sim_functions_pedigreed.R")

parpar = data.frame(n.pop0 = 100,
                    end.time = 15,
                    init.row = 1e4,
                    n.loci = 25,
                    w.max = 2,
                    theta = 3,
                    wfitn = sqrt(1 / 0.14 / 2),
                    sig.e = 0.5,
                    pos.p = 0.5,
                    alpha = 0.000)

set.seed(4199)

testsim = sim.par.id(a = c(1/2, -1/2), params = parpar)

table(testsim$gen)

zh_t = 2:parpar$end.time

for (timestep in 1:(parpar$end.time-1)) {
  gens = testsim %>%
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
  
  #gh_t[timestep] = summary(lm(g_i ~ g_midp, midpars))$coefficients[,'Estimate']['g_midp']
  zh_t[timestep] = summary(lm(z_i ~ z_midp, midpars))$coefficients[,'Estimate']['z_midp']
}

plot(zh_t, type = 'l')
mean(zh_t)

# this is not that far from a predicted z of 0.67...

actual.v = testsim %>%
  group_by(gen) %>%
  summarise(v = var(g_i))

predic.v = testsim %>%
  select(-c(i, g_i, z_i, w_i, r_i, fem, mom_i, dad_i)) %>%
  gather(key = loc.copy, value = val, -c(gen)) %>%
  mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
  group_by(gen, locus) %>%
  summarise(p = mean(val > 0)) %>%
  group_by(gen) %>%
  summarise(v = sum(2 * p * (1 - p)) / parpar$n.loci[1]) 

plot(predic.v$v, actual.v$v)
abline(0,1)

# The estimate for genetic variation looks right to me.

testsim %>%
  mutate(E = g_i - z_i) %>%
  ggplot() +
  geom_density(aes(x = E, group = gen, colour = factor(gen))) +
  scale_color_viridis_d()


testsim %>%
  group_by(gen) %>%
  mutate(d_i = parpar$theta - z_i) %>%
  summarise(dbar = mean(d_i),
            # sig2 = var(g_i),
            dwbar = mean(d_i * w_i / mean(w_i)),
            preddw = dbar * with(parpar, wfitn^2 / (wfitn^2 + var(z_i))),
            h_t = var(g_i) / var(z_i),
            preddbar1 = dbar * with(parpar, (wfitn^2 + sig.e^2) / (wfitn^2 + sig.e^2 + var(g_i))),
            preddbar2 = h_t * dwbar + (1-h_t) * dbar,
            preddbar3 = (1/3) * dwbar + (2/3) * dbar) %>%
  mutate(dbar_tp1 = c(dbar[-1], 0)) %>%
  filter(gen < 15) %>%
  ggplot() +
  geom_point(aes(x = preddbar3, y = dbar_tp1)) +
  geom_segment(aes(x = min(preddbar3), xend = max(preddbar3),
                   y = min(preddbar3), yend = max(preddbar3)))

# one thing that seems amiss... the h_t is kind of low
# but putting in the approximate correct heritability (preddbar3) doesn't remedy things
# (but try swapping h_t and 1-h_t and the estimates look better...?
# this actually looks right...?
# is this a falsity??



### aghghh what the heck is happening???
# (see above)

# Phenotypes appear to be correctly coming from genotypes
# The lande expression for zwt works (fitness-weighted phenotypes work)
# Heritability, I think, is working correctly...
