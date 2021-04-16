library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)

all.data = 
  rbind(
    read.csv('simulations/outputs/longsims/longsims_n100_a000_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n100_a000_lowvar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a000_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a000_lowvar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n100_a035_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n100_a035_lowvar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a035_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a035_lowvar.csv')
  )

# write.csv(all.data, row.names = FALSE,
#           'simulations/outputs/final_results/longsims/longsims_combined.csv')

nrow(all.data)

hmm = all.data %>%
  arrange(trial, n.pop0, low.var, alpha, desc(gen)) %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  filter(gen < 50)

head(hmm)

hmm %>% ggplot(aes(x = n)) + geom_histogram() + facet_wrap(n.pop0 + low.var ~ alpha, nrow = 2)

all.data = all.data %>%
  arrange(trial, n.pop0, low.var, alpha, desc(gen)) %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  group_by(n.pop0, low.var, alpha, trial) %>%
  summarise(extinct = gen < 50 & n < 10000) %>%
  merge(x = all.data, y = ., by = c('trial', 'n.pop0', 'low.var', 'alpha'))

# Now do some stuff

# Extinction probability:

all.data %>%
  distinct(n.pop0, low.var, alpha, trial, .keep_all = TRUE) %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(n = mean(extinct))

all.ext = all.data %>%
  group_by(n.pop0, low.var, alpha, trial) %>%
  summarise(ext.gen = ifelse(extinct[1], max(gen), 50))

all.ext %>%
  filter(ext.gen < 50) %>%
  ungroup() %>%
  mutate(n.pop0 = factor(ifelse(n.pop0 %in% 100, "Large", "Small")),
         alpha = factor(ifelse(alpha > 0, "Density dependent", "Density independent")),
         low.var = factor(ifelse(low.var, "Low diversity", "High diversity"))) %>%
  ggplot() +
  geom_histogram(aes(x = ext.gen, fill = factor(alpha)), 
                 position = 'identity',
                 binwidth = 1, alpha = 0.5) +
  scale_fill_manual(values = c('purple', 'black')) +
  facet_wrap(~ paste(n.pop0, low.var, sep = ', '), nrow = 4) +
  theme(legend.position = 'none')

all.ext %>%
  filter(ext.gen < 50) %>%
  ggplot() +
  geom_density(aes(x = ext.gen, colour = factor(alpha))) +
  scale_colour_manual(values = c('black', 'purple')) +
  facet_wrap(n.pop0 ~ low.var) +
  theme(legend.position = 'none')

all.ext %>%
  filter(ext.gen < 50) %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(tau = mean(ext.gen),
            tau.se = sd(ext.gen / n())) # not quite the way to get s.e. for this but sufficient

all.ext %>%
  filter(ext.gen < 50) %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(tau = mean(ext.gen)) %>%
  spread(key = alpha, value = tau)

all.ext %>% filter(ext.gen < 50) %>%
  glm(ext.gen ~ factor(n.pop0) * low.var * factor(alpha), data = ., family = 'poisson') %>%
  summary()

all.ext %>% filter(ext.gen < 50) %>%
  glm(ext.gen ~ factor(n.pop0) + low.var + 
                factor(n.pop0) * factor(alpha) + 
                low.var * factor(alpha), 
      data = ., family = 'poisson') %>%
  summary()

all.ext %>% filter(ext.gen < 50) %>%
  glm(ext.gen ~ factor(n.pop0) + low.var + factor(n.pop0) * factor(alpha), 
      data = ., family = 'poisson') %>%
  summary()
  
# Feel like a censored approach here would be helpful!

### Plot accumulation of extinctions
all.extinctions.long = all.data %>%
  select(gen, trial, n.pop0, low.var, alpha, extinct) %>%
  group_by(trial, n.pop0, low.var, alpha) %>%
  mutate(ext.this.gen = (gen %in% max(gen)) & extinct) %>%
  group_by(n.pop0, low.var, alpha, gen) %>%
  summarise(n.extinctions = sum(ext.this.gen)) %>%
  ungroup() %>%
  arrange(n.pop0, low.var, alpha, gen) %>%
  mutate(n.pop0 = factor(ifelse(n.pop0 %in% 100, "Large", "Small")),
         alpha = factor(ifelse(alpha > 0, "Density dependent", "Density independent")),
         low.var = factor(ifelse(low.var, "Low diversity", "High diversity")))

cuml.plot = all.extinctions.long %>%
  group_by(n.pop0, low.var, alpha) %>%
  mutate(p.extinct = cumsum(n.extinctions / 1000)) %>%
  ggplot(aes(x = gen)) +
  geom_ribbon(
    aes(
      xmin = 0, xmax = 50,
      ymin = 0, ymax = p.extinct,
      group = interaction(n.pop0, low.var, alpha),
      fill = factor(alpha)
    ),
    alpha = 0.5
  ) +
  scale_fill_manual(values = c('purple', 'black')) +
  labs(x = 'Generation', y = 'Probability of extinction') +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        strip.background = element_rect(colour = 'black'))

inst.plot = all.extinctions.long %>%
  group_by(n.pop0, low.var, alpha) %>%
  mutate(n.extant = 1000 - c(0, cumsum(n.extinctions)[-50]),
         p.instant.extinct = n.extinctions / n.extant) %>%
  ggplot(aes(x = gen)) +
  geom_ribbon(
    aes(
      ymin = p.instant.extinct - 2 * sqrt(p.instant.extinct*(1-p.instant.extinct) / n.extant),
      ymax = p.instant.extinct + 2 * sqrt(p.instant.extinct*(1-p.instant.extinct) / n.extant),
      fill = factor(alpha)
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      y = p.instant.extinct,
      group = factor(alpha),
      colour = factor(alpha)
    ),
    size = 1
  ) +
  labs(x = 'Generation', y = '') +
  scale_fill_manual(values = c('purple', 'black')) +
  scale_color_manual(values = c('purple', 'black')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        strip.background = element_rect(colour = 'black'))

all.data %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  glm(extinct ~ factor(n.pop0) * low.var * factor(alpha) * gbar,
      data = ., family = 'binomial') %>%
  summary()

all.data %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  glm(extinct ~ factor(n.pop0) * low.var * factor(alpha) +
                factor(n.pop0) * low.var * gbar +
                factor(n.pop0) * gbar * factor(alpha) +
                gbar * low.var * factor(alpha),
      data = ., family = 'binomial') %>%
  summary()

all.data %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  glm(extinct ~ factor(n.pop0) * low.var * factor(alpha) +
                factor(n.pop0) * low.var * gbar +
                factor(n.pop0) * gbar * factor(alpha),
      data = ., family = 'binomial') %>%
  summary()

all.data %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  glm(extinct ~ factor(n.pop0) * low.var * factor(alpha) +
                factor(n.pop0) * gbar * factor(alpha) +
                low.var * gbar,
      data = ., family = 'binomial') %>%
  summary()

all.data %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  glm(extinct ~ factor(n.pop0) * low.var * factor(alpha) +
                factor(n.pop0) * gbar +
                gbar * factor(alpha) +
                low.var * gbar,
      data = ., family = 'binomial') %>%
  summary()

all.data %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  glm(extinct ~ factor(n.pop0) * low.var * factor(alpha) +
                gbar * factor(alpha) +
                low.var * gbar,
      data = ., family = 'binomial') %>%
  summary()

all.data %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  glm(extinct ~ factor(n.pop0) * low.var * factor(alpha) +
                low.var * gbar,
      data = ., family = 'binomial') %>%
  summary()

# kinda puny model!

ext.mod.final =  glm(
  extinct ~ factor(n.pop0) * low.var * factor(alpha) +
            low.var * gbar,
  data = all.data %>% distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE), 
  family = 'binomial'
  ) 

geno.ext.preds = expand.grid(gbar = (-5:5)/10,
            alpha = c(0, 0.0035),
            n.pop0 = factor(c(20, 100)),
            low.var = c(TRUE, FALSE)) %>%
  (function(x) cbind(x, p = predict(ext.mod.final, x, type = 'response'))) %>%
  mutate(n.pop0 = factor(ifelse(n.pop0 %in% 100, "Large", "Small")),
         alpha = factor(ifelse(alpha > 0, "Density dependent", "Density independent")),
         low.var = factor(ifelse(low.var, "Low diversity", "High diversity")))

geno.plot = geno.ext.preds %>%
  ggplot() +
  geom_line(
    aes(
      x = gbar, y = p, group = factor(alpha),
      colour = factor(alpha)
    )
  ) +
  # geom_point(
  #   data = all.data %>% 
  #     distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  #     mutate(ext.jitter = as.numeric(extinct) + (1 - 2*as.numeric(extinct)) * runif(nrow(.), 0, 0.1)),
  #   aes(
  #     x = gbar, y = ext.jitter, colour = factor(alpha)
  #   ),
  #   alpha = 0.1
  # ) +
  labs(x = 'Initial genotype', y = '') +
  scale_color_manual(values = c('purple', 'black')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        strip.background = element_rect(colour = 'black'))


extinct.legend = get_legend(
  inst.plot + 
    guides(fill = guide_legend(''), colour = guide_legend('')) +
    theme(legend.position = 'bottom',
          legend.box.margin = margin(1, 0, 0, 0))
)

data.plots = plot_grid(cuml.plot, inst.plot, geno.plot, 
                       labels = c('(A)', '(B)', '(C)'),
                       nrow = 1)
plot_grid(data.plots, extinct.legend, ncol = 1, rel_heights = c(1, .1)) %>%
  save_plot(filename = 'simulations/analysis_results/figure_drafts/draft_figs/fig_2_long.pdf',
            base_width = 8, base_height = 8)

# Plotting all time serieses 

all.data %>%
  ggplot() +
  geom_line(aes(x = gen, y = n, colour = extinct, group = trial),
            size = 0.2) +
  scale_y_log10() +
  scale_color_manual(values = c('black', 'red')) +
  facet_wrap(n.pop0 + low.var ~ alpha, ncol = 4)

all.data %>%
  filter(trial < 26) %>%
  ggplot() +
  geom_line(aes(x = gen, y = n, colour = factor(alpha), 
                linetype = extinct, group = interaction(alpha, trial)),
            size = 0.4) +
  scale_y_log10() +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap(n.pop0 + low.var ~ ., ncol = 1) +
  theme(legend.position = 'none')

all.data %>%
  filter(trial < 26) %>%
  mutate(n.pop0 = factor(ifelse(n.pop0 %in% 100, "Large", "Small")),
         alpha = factor(ifelse(alpha > 0, "Density dependent", "Density independent")),
         low.var = factor(ifelse(low.var, "Low diversity", "High diversity"))) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      x = gen, 
      y = n, 
      colour = factor(alpha), 
      linetype = extinct, 
      group = interaction(alpha, trial)
    ),
  size = 0.6) +
  geom_line(
    aes(
      x = gen, 
      y = n, 
      colour = factor(alpha), 
      linetype = extinct, 
      group = interaction(alpha, trial),
      alpha = factor(alpha)
    ),
  size = 0.6) +
  labs(x = 'Generation', y = '') +
  scale_y_log10() +
  scale_alpha_manual(values = c(0, 1)) +
  scale_fill_manual(values = c('purple', 'black')) +
  scale_color_manual(values = c('purple', 'black')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        strip.background = element_rect(colour = 'black'))

# Do some visualizations of g, v for extinction risk

all.data %>%
  ggplot() +
  geom_line(aes(x = gen, y = gbar, 
                colour = extinct, group = trial),
            size = 0.2) +
  scale_color_manual(values = c('black', 'red')) +
  facet_wrap(n.pop0 + low.var ~ alpha, ncol = 4) +
  theme(legend.position = 'none')

all.data %>%
  filter(trial < 200) %>%
  ggplot() +
  geom_line(aes(x = gen, y = wbar, 
                colour = extinct, group = trial),
            size = 0.2) +
  scale_color_manual(values = c('black', 'red')) +
  facet_wrap(n.pop0 + low.var ~ alpha, ncol = 4) +
  theme(legend.position = 'none')
  
all.data %>%
  filter(trial < 200) %>%
  ggplot() +
  geom_line(aes(x = gen, y = v, 
                colour = extinct, group = trial),
            size = 0.2) +
  scale_color_manual(values = c('black', 'red')) +
  facet_wrap(n.pop0 + low.var ~ alpha, ncol = 4) +
  theme(legend.position = 'none')

# Okay... FUNctionl programming time folks

library(purrr)

all.tau = all.data %>%
  group_by(trial, n.pop0, low.var, alpha) %>%
  mutate(loglam = c(NA, diff(log(n))),
         d = abs(2.75 - gbar),
         k = c(NA, exp(diff(log(d)))),
         varred = c(NA, exp(diff(log(v))))) %>%
  mutate(tau = ifelse(extinct, max(gen) - gen + 1 - as.numeric(any(n < 2)), 50),
         ext = extinct | any(n < 2),
         ext5 = tau < 5) %>%
  filter(!is.na(loglam), gen > 1) %>%
  ungroup() %>%
  mutate(n.pop0 = factor(n.pop0),
         alpha  = factor(alpha))

ll.blank = expand.grid(
  loglam = (-7:4) / 10,
  low.var = c(TRUE, FALSE),
  n.pop0 = factor(c(20, 100)),
  alpha = factor(c(0, 0.0035))
)

llmods = all.tau %>%
  filter(gen < 25) %>%
  split(.$gen) %>%
  map(~ glm(ext5 ~ low.var * n.pop0 * alpha * loglam,
            data = .,
            family = 'binomial'))

ll.preds = lapply(llmods, function(x) predict.glm(object = x, newdata = ll.blank, type = 'resp')) %>%
  do.call(what = cbind) %>%
  cbind(ll.blank, .) %>%
  gather(key = gen, value = p.ext, -c(loglam, low.var, n.pop0, alpha)) %>%
  mutate(gen = as.numeric(gen))

ggplot(ll.preds) +
  geom_tile(
    aes(
      x = gen,
      y = loglam,
      fill = p.ext
    )
  ) +
  geom_segment(
    aes(x = 1, xend = 25, y = 0, yend = 0),
    colour = 'gray',
    size = 0.2
  ) +
  facet_wrap(paste(n.pop0, low.var) ~ alpha, nrow = 4) +
  scale_fill_viridis_b(option = 'B') +
  theme(legend.position = 'bottom',
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'black'))

###

# Assessing time until original size is reached (>3 generations)

which.diff = function(x, lag) {
  if(length(x) >= lag) { x[which.min(diff(x, lag = (lag - 1)) == (lag-1))] 
  } else { 0 }
}

rescue.summ = all.data %>%
  filter(n > n.pop0, gen > 1) %>%
  group_by(trial, n.pop0, low.var, alpha) %>%
  summarise(rescue.time = which.diff(gen, 3)) %>%
  filter(rescue.time > 0)

rescue.summ %>%
  ggplot() +
  geom_histogram(aes(x = rescue.time, fill = factor(alpha)), 
                 position = 'identity', alpha = 0.5,
                 binwidth = 1) +
  facet_wrap(n.pop0 ~ low.var) 

rescue.summ %>%
  ggplot() +
  geom_histogram(aes(x = rescue.time, fill = factor(alpha)), 
                 position = 'identity', alpha = 0.5,
                 binwidth = 0.1) +
  scale_x_log10() +
  facet_wrap(n.pop0 ~ low.var) 

# rescue probably more likely in small populations with a large amount of drift
# confine to the middle 80% of the distributon of large population genotypes

rescue.summ.censor = all.data %>%
  # Getting the middle 80th percentiles of the large populations in each group
  filter(n.pop0 > 20) %>%
  distinct(trial, low.var, alpha, .keep_all = TRUE) %>%
  group_by(low.var, alpha) %>%
  summarise(g10 = quantile(gbar, 0.1),
            g90 = quantile(gbar, 0.9)) %>%
  ungroup() %>%
  merge(y = all.data,
        by = c("low.var", "alpha")) %>%
  # Filter out the outlier populations
  group_by(trial, n.pop0, low.var, alpha) %>%
  filter(gbar[1] < g90 & gbar[1] > g10) %>%
  select(-c(g10, g90)) %>%
  # Estimate the rescue times
  filter(n > n.pop0, gen > 1) %>%
  group_by(trial, n.pop0, low.var, alpha) %>%
  summarise(rescue.time = which.diff(gen, 3)) %>%
  filter(rescue.time > 0)

nrow(rescue.summ) - nrow(rescue.summ.censor)

rescue.summ.censor %>%
  ggplot() +
  geom_histogram(aes(x = rescue.time, fill = factor(alpha)), 
                 position = 'identity', alpha = 0.5,
                 binwidth = 1) +
  facet_wrap(n.pop0 ~ low.var) 

rescue.summ.censor %>%
  ggplot() +
  geom_histogram(aes(x = rescue.time, fill = factor(alpha)), 
                 position = 'identity', alpha = 0.5,
                 binwidth = 0.1) +
  scale_x_log10() +
  facet_wrap(n.pop0 ~ low.var) 

merge(x = all.data %>% distinct(trial, n.pop0, low.var, alpha),
      y = rescue.summ %>% ungroup(),
      by = c('trial', 'n.pop0', 'low.var', 'alpha'),
      all.x = TRUE, all.y = TRUE) %>%
  mutate(rescued = !is.na(rescue.time)) %>%
  glm(formula = rescued ~ factor(n.pop0) * low.var * factor(alpha) - 1,
      family = 'quasibinomial') %>%
  summary()

wescue.summ = all.data %>%
  mutate(rbar = wbar * exp(-alpha * n)) %>%
  filter(rbar > 1, gen > 1) %>%
  group_by(trial, n.pop0, low.var, alpha) %>%
  summarise(rescue.time = which.diff(gen, 3)) %>%
  filter(rescue.time > 0)

wescue.summ

wescue.summ %>%
  ggplot() +
  geom_histogram(aes(x = rescue.time, fill = factor(alpha)), 
                 position = 'identity', alpha = 0.5,
                 binwidth = 1) +
  facet_wrap(n.pop0 ~ low.var) 


merge(x = all.data %>% distinct(trial, n.pop0, low.var, alpha),
      y = wescue.summ %>% ungroup(),
      by = c('trial', 'n.pop0', 'low.var', 'alpha'),
      all.x = TRUE, all.y = TRUE) %>%
  mutate(rescued = !is.na(rescue.time)) %>%
  glm(formula = rescued ~ factor(n.pop0) * low.var * factor(alpha) - 1,
      family = 'quasibinomial') %>%
  summary()

wescue.summ %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(p = n() / 1000,
            g = mean(rescue.time))

### Rescue mods with initial genotypes

p.wescue = merge(x = all.data %>% 
                   distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
                   select(trial, n.pop0, low.var, alpha, gbar),
                 y = wescue.summ %>% ungroup(),
                 by = c('trial', 'n.pop0', 'low.var', 'alpha'),
                 all.x = TRUE, all.y = TRUE) %>%
  mutate(rescued = !is.na(rescue.time))

# Four way model
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) * gbar,
    data = p.wescue,
    family = 'quasibinomial') %>%
  summary()

# All three-way combos
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) +
              factor(n.pop0) * low.var * gbar +
              factor(n.pop0) * factor(alpha) * gbar +
              low.var * factor(alpha) * gbar,
    data = p.wescue,
    family = 'quasibinomial') %>%
  summary()

# Remove variance-alpha-gbar term
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) +
              factor(n.pop0) * low.var * gbar +
              factor(n.pop0) * factor(alpha) * gbar,
    data = p.wescue,
    family = 'quasibinomial') %>%
  summary()

# Remove size-variance-gbar
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) +
              factor(n.pop0) * factor(alpha) * gbar +
              low.var * gbar,
    data = p.wescue,
    family = 'quasibinomial') %>%
  summary()

# Remove size-alpha-gbar
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) +
              factor(n.pop0) * gbar +
              factor(alpha) * gbar +
              low.var * gbar,
    data = p.wescue,
    family = 'quasibinomial') %>%
  summary()

# Remove size-gbar 
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) +
              factor(alpha) * gbar +
              low.var * gbar,
    data = p.wescue,
    family = 'quasibinomial') %>%
  summary()

# Remove alpha-gbar
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) +
              low.var * gbar,
    data = p.wescue,
    family = 'quasibinomial') %>%
  summary()

# At alpha of 0.1, remove lowvar-gbar
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) + gbar,
    data = p.wescue,
    family = 'quasibinomial') %>%
  summary()
# This is good.

final.wescue.mod = glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) + gbar,
                       data = p.wescue,
                       family = 'quasibinomial')

# Wescue genotype plots:

geno.wes.preds = expand.grid(gbar = (-5:5)/10,
                             alpha = c(0, 0.0035),
                             n.pop0 = factor(c(20, 100)),
                             low.var = c(TRUE, FALSE)) %>%
  (function(x) cbind(x, p = predict(final.wescue.mod, x, type = 'response'))) %>%
  mutate(n.pop0 = factor(ifelse(n.pop0 %in% 100, "Large", "Small")),
         alpha = factor(ifelse(alpha > 0, "Density dependent", "Density independent")),
         low.var = factor(ifelse(low.var, "Low diversity", "High diversity")))

geno.wes.preds %>%
  ggplot() +
  geom_line(
    aes(
      x = gbar, y = p, group = factor(alpha),
      colour = factor(alpha)
    )
  ) +
  labs(x = 'Initial genotype', y = '') +
  scale_color_manual(values = c('purple', 'black')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        strip.background = element_rect(colour = 'black'))

### Do same plot but for rescue (original size reached)
# (results should be very similar)

p.rescue = merge(x = all.data %>% 
                   distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
                   select(trial, n.pop0, low.var, alpha, gbar),
                 y = rescue.summ %>% ungroup(),
                 by = c('trial', 'n.pop0', 'low.var', 'alpha'),
                 all.x = TRUE, all.y = TRUE) %>%
  mutate(rescued = !is.na(rescue.time))


# Four way model
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) * gbar,
    data = p.rescue,
    family = 'binomial') %>%
  summary()

# All three-way combos
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) +
              factor(n.pop0) * low.var * gbar +
              factor(n.pop0) * factor(alpha) * gbar +
              low.var * factor(alpha) * gbar,
    data = p.rescue,
    family = 'binomial') %>%
  summary()

# All three-way combos
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) +
              factor(n.pop0) * low.var * gbar +
              factor(n.pop0) * factor(alpha) * gbar +
      low.var * factor(alpha) * gbar,
    data = p.rescue,
    family = 'binomial') %>%
  summary()

# Remove var-alpha-gbar
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) +
              factor(n.pop0) * low.var * gbar +
              factor(n.pop0) * factor(alpha) * gbar,
    data = p.rescue,
    family = 'binomial') %>%
  summary()

# Remove size-var-gbar (p = 0.035, alpha = 0.01)
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) +
              factor(n.pop0) * factor(alpha) * gbar +
              low.var * gbar,
    data = p.rescue,
    family = 'binomial') %>%
  summary()

# Remove size-alpha-gbar (p = 0.030, alpha = 0.01)
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) +
              factor(n.pop0) * gbar +
              factor(alpha) * gbar +
              low.var * gbar,
    data = p.rescue,
    family = 'binomial') %>%
  summary()

# Remove size-gbar
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) +
              factor(alpha) * gbar +
              low.var * gbar,
    data = p.rescue,
    family = 'binomial') %>%
  summary()

# Remove alpha-gbar
glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) +
              low.var * gbar,
    data = p.rescue,
    family = 'binomial') %>%
  summary()

# This is it, chief.
final.rescue.mod = glm(rescued ~ factor(n.pop0) * low.var * factor(alpha) +
                                 low.var * gbar,
                       data = p.rescue,
                       family = 'binomial')

# Make the predictions
geno.res.preds = expand.grid(gbar = (-5:5)/10,
                             alpha = c(0, 0.0035),
                             n.pop0 = factor(c(20, 100)),
                             low.var = c(TRUE, FALSE)) %>%
  (function(x) cbind(x, p = predict(final.rescue.mod, x, type = 'response'))) %>%
  mutate(n.pop0 = factor(ifelse(n.pop0 %in% 100, "Large", "Small")),
         alpha = factor(ifelse(alpha > 0, "Density dependent", "Density independent")),
         low.var = factor(ifelse(low.var, "Low diversity", "High diversity")))

# Make the plot.
geno.res.preds %>%
  ggplot() +
  geom_line(
    aes(
      x = gbar, y = p, group = factor(alpha),
      colour = factor(alpha)
    )
  ) +
  labs(x = 'Initial genotype', y = '') +
  scale_color_manual(values = c('purple', 'black')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        strip.background = element_rect(colour = 'black'))

# Combined:
merge(geno.res.preds, geno.wes.preds,
      by = c('gbar', 'alpha', 'n.pop0', 'low.var'),
      suffixes = c('.size', '.fitn')) %>%
  ggplot() +
  geom_line(
    aes(
      x = gbar, y = p.fitn, group = factor(alpha),
      colour = factor(alpha)
    ),
    linetype = 2, size = 1.1
  ) +
  geom_line(
    aes(
      x = gbar, y = p.size, group = factor(alpha),
      colour = factor(alpha)
    ),
    linetype = 1
  ) +
  labs(x = 'Initial genotype', y = '') +
  scale_color_manual(values = c('purple', 'black')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 2) +
  theme(legend.position = 'none',
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        strip.background = element_rect(colour = 'black'))
