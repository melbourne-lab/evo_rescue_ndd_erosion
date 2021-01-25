library(ggplot2)

rm(list = ls())

source('base_model_source/sim_functions.R')
source('base_model_source/sim_functions_overlap.R')

pars = data.frame(
  end.time = 10,
  init.row = 1e4,
  n.loci = 25,
  n.pop0 = 50,
  w.max = 2,
  theta = 3,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0
)

n.sims = 400

list.nonover = vector(mode = 'list', length = n.sims)
list.overlap = vector(mode = 'list', length = n.sims)

set.seed(29)

for (trial in 1:n.sims) {
  
  list.nonover[[trial]] = sim(params = pars)
  
  list.overlap[[trial]] = sim.overlap(params = pars)
  
  print(trial)
  
}

# Note: this is super slow! 

sims.no = unroller(list.nonover)
sims.ov = unroller(list.overlap)

# sims.no %>%
#   group_by(trial, gen) %>%
#   summarise(n_kt = n()) %>%
#   group_by(gen) %>%
#   summarise(n_t = mean(n_kt)) %>%
#   ggplot() +
#   geom_line(aes(x = gen, y = n_t))
# 
# 
# sims.ov %>%
#   group_by(trial, gen) %>%
#   summarise(n_kt = n()) %>%
#   group_by(gen) %>%
#   summarise(n_t = mean(n_kt)) %>%
#   ggplot() +
#   geom_line(aes(x = gen, y = n_t))

# Compile into one data frame

sims.all = rbind(
  sims.no %>% mutate(overlap = FALSE),
  sims.ov %>% mutate(overlap = TRUE) %>% select(-age)
)

write.csv(sims.all, file = 'analyze_model/overlapping_gens/sims_800.csv',
          row.names = FALSE)

# Compare extinction rates
sims.all %>%
  distinct(trial, gen, overlap) %>%
  group_by(gen, overlap) %>%
  summarise(p.extant = n() / n.sims) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = p.extant, group = overlap, colour = overlap))

# Compare population sizes
# (note: this is extant only, subject to extnction bias)
sims.all %>%
  group_by(overlap, gen, trial) %>%
  summarise(n_t = n()) %>%
  ungroup() %>%
  rbind(
    expand.grid(
      overlap = c(TRUE, FALSE),
      gen = 1:pars$end.time,
      trial = 1:n.sims,
      n_t = 0
    )
  ) %>%
  group_by(overlap, gen, trial) %>%
  summarise(n_t = sum(n_t)) %>%
  group_by(overlap, gen) %>%
  summarise(nbar = mean(n_t),
            nvar = var(n_t),
            nobs = n()) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = nbar,
      group = overlap,
      color = overlap
    )
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / nobs),
      ymax = nbar + 2 * sqrt(nvar / nobs),
      group = overlap,
      fill = overlap
    ),
    alpha = 0.5
  )

sims.all %>%
  group_by(overlap, gen, trial) %>%
  summarise(g_t = mean(g_i)) %>%
  group_by(overlap, gen) %>%
  summarise(gbar = mean(g_t),
            gvar = var(g_t),
            nobs = n()) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = gbar,
      group = overlap,
      color = overlap
    )
  ) +
  geom_ribbon(
    aes(
      ymin = gbar - 2 * sqrt(gvar / nobs),
      ymax = gbar + 2 * sqrt(gvar / nobs),
      group = overlap,
      fill = overlap
    ),
    alpha = 0.5
  )

