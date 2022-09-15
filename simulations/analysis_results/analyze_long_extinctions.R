### Script containing all results pertaining to extinctions in the long (50
### generation) simulations.
### Contains code for Figure 2
### SN - 29 Mar 2021


# Load packages
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)

### Setup and declaring variables

# Clear namespace
rm(list = ls())

# Load in simulation data
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

# Add in extinction designation
all.data = all.data %>%
  arrange(trial, n.pop0, low.var, alpha, desc(gen)) %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  group_by(n.pop0, low.var, alpha, trial) %>%
  summarise(extinct = gen < 50 & n < 10000) %>%
  merge(x = all.data, y = ., by = c('trial', 'n.pop0', 'low.var', 'alpha'))

# Get number of extinctions occurring in each generation
all.extinctions.long = all.data %>%
  select(gen, trial, n.pop0, low.var, alpha, extinct) %>%
  group_by(trial, n.pop0, low.var, alpha) %>%
  mutate(ext.this.gen = (gen %in% max(gen)) & extinct) %>%
  group_by(n.pop0, low.var, alpha, gen) %>%
  summarise(n.extinctions = sum(ext.this.gen)) %>%
  ungroup() %>%
  arrange(n.pop0, low.var, alpha, gen) %>%
  mutate(romnum = factor(ifelse(n.pop0 %in% 100,
                                ifelse(low.var, '(i)', '(ii)'),
                                ifelse(low.var, '(iii)', '(iv)'))),
         n.pop0 = factor(ifelse(n.pop0 %in% 100, "Large", "Small")),
         alpha = factor(ifelse(alpha > 0, "Density dependent", "Density independent")),
         low.var = factor(ifelse(low.var, "Low diversity", "High diversity")))

# Run extinction model
ext.data = all.data %>% distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE)

ext.mod.final =  glm(
  extinct ~ factor(n.pop0) * low.var * factor(alpha) +
    low.var * gbar,
  data = ext.data, 
  family = 'binomial'
)

ext.mod.full =  glm(
  extinct ~ factor(n.pop0) * low.var * factor(alpha) * gbar,
  data = ext.data, 
  family = 'binomial'
)

summary(ext.mod.final)
summary(ext.mod.full)

### Make plots

# # Make object for labelling plots
# pan.labls = expand.grid(low.var = factor(c('High diversity', 'Low diversity')),
#             n.pop0 = factor(c('Large', 'Small'))) %>%
#   mutate(labl = c('(i)', '(ii)', '(iii)', '(iv)'))

# Instantaneous probability of extinction (a)
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
  # geom_text(
  #   data = pan.labls,
  #   aes(
  #     x = 1, y = 0.6,
  #     label = labl
  #   )
  # ) +
  labs(x = 'Generation', y = '') +
  scale_fill_manual(values = c('purple', 'black')) +
  scale_color_manual(values = c('purple', 'black')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        strip.background = element_rect(colour = 'black'))

# Cumulative extinction plot (b)
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
  # geom_text(
  #   data = pan.labls,
  #   aes(
  #     x = 1, y = 0.7,
  #     label = labl
  #   )
  # ) +
  scale_fill_manual(values = c('purple', 'black')) +
  labs(x = 'Generation', y = 'Probability of extinction') +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        strip.background = element_rect(colour = 'black'))

# Make predicted extinctions based on genotype
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
  # geom_text(
  #   data = pan.labls,
  #   aes(
  #     x = 0, y = 0.7,
  #     label = labl
  #   )
  # ) +
  labs(x = 'Initial genotype', y = '') +
  scale_color_manual(values = c('purple', 'black')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        strip.background = element_rect(colour = 'black'))

geno.ext.fulls = expand.grid(gbar = (-5:5)/10,
                             alpha = c(0, 0.0035),
                             n.pop0 = factor(c(20, 100)),
                             low.var = c(TRUE, FALSE)) %>%
  (function(x) cbind(x, p = predict(ext.mod.full, x, type = 'response'))) %>%
  mutate(n.pop0 = factor(ifelse(n.pop0 %in% 100, "Large", "Small")),
         alpha = factor(ifelse(alpha > 0, "Density dependent", "Density independent")),
         low.var = factor(ifelse(low.var, "Low diversity", "High diversity")))

geno.full = geno.ext.fulls %>%
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
  # geom_text(
  #   data = pan.labls,
#   aes(
#     x = 0, y = 0.7,
#     label = labl
#   )
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

data.plots = plot_grid(cuml.plot, inst.plot, geno.full, 
                       labels = c('(A)', '(B)', '(C)'),
                       nrow = 1)

plot_grid(data.plots, extinct.legend, ncol = 1, rel_heights = c(1, .1)) %>%
  save_plot(filename = 'simulations/analysis_results/figure_drafts/draft_figs/fig_2_long.pdf',
            base_width = 8, base_height = 8)

### Some summary statistics

# Increased extinction risk for each parameter combo
ext.data %>% 
  group_by(n.pop0, low.var, alpha) %>%
  summarise(pext = mean(extinct)) %>%
  spread(key = alpha, value = pext) %>%
  mutate(incr = `0.0035` / `0`)

### Try a figure with opposite direction of paneling

size.plot.row = all.data %>%
  filter(trial < 26) %>%
  mutate(romnum = factor(ifelse(n.pop0 %in% 100,
                                ifelse(!low.var, '(i)', '(ii)'),
                                ifelse(!low.var, '(iii)', '(iv)'))),
         n.pop0 = factor(ifelse(n.pop0 %in% 100, "Large", "Small")),
         alpha = factor(ifelse(alpha > 0, "Density dependent", "Density independent")),
         low.var = factor(ifelse(!low.var, "High diversity", "Low diversity"))) %>%
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
    size = 0.5) +
  labs(x = 'Generation', y = '') +
  scale_y_log10() +
  scale_alpha_manual(values = c(0, 1)) +
  scale_fill_manual(values = c('purple', 'black')) +
  scale_color_manual(values = c('purple', 'black')) +
  facet_wrap( ~ paste0(romnum, ' ', n.pop0, ', ', low.var), ncol = 4) +
  theme(legend.position = 'none',
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        strip.background = element_rect(colour = 'black'),
        strip.text = element_text(size = 11))

inst.plot.row = all.extinctions.long %>%
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
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        plot.margin = margin(b = 0, r = 5, l = 20, unit = 'pt'))

# Cumulative extinction plot (b)
cuml.plot.row = all.extinctions.long %>%
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
  labs(x = 'Generation', y = '') +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  theme(legend.position = 'none',        
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        plot.margin = margin(b = 0, r = 5, l = 20, unit = 'pt'))

# geno.full.row = geno.ext.fulls %>%
#   ggplot() +
#   geom_line(
#     aes(
#       x = gbar, y = p, group = factor(alpha),
#       colour = factor(alpha)
#     )
#   ) +
#   labs(x = 'Initial genotype', y = '') +
#   scale_color_manual(values = c('purple', 'black')) +
#   facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
#   theme(legend.position = 'none',
#         strip.text = element_blank(),
#         strip.background = element_blank(),
#         panel.grid.major = element_line(colour = 'gray88'),
#         panel.background = element_rect(fill = 'white'),
#         plot.margin = margin(b = 0, r = 5, l = 5, unit = 'pt'))

data.plots.row = plot_grid(size.plot.row, cuml.plot.row, 
                           inst.plot.row, geno.full.row, 
                           labels = c('(A)', '(B)', '(C)', '(D)'),
                           nrow = 4)

plot_grid(data.plots.row, extinct.legend, ncol = 1, rel_heights = c(1, .1)) %>%
  save_plot(filename = 'simulations/analysis_results/figure_drafts/draft_figs/fig_2_long_rows.pdf',
            base_width = 8, base_height = 8)

##### Bayesian model

library(rstanarm)

ext.data = ext.data %>%
  mutate(low.var = factor(ifelse(low.var, 'low.div', 'hi.div')),
         n.pop0  = factor(ifelse(n.pop0 > 50, 'large', 'small')),
         alpha   = factor(ifelse(alpha > 0, 'ndd', 'di')))

ext.mod.full =  stan_glm(
  extinct ~ n.pop0 * low.var * alpha * gbar,
  data = ext.data, 
  seed = 992,
  cores = 4,
  prior = normal(0, 5),
  family = 'binomial'
)

summary(ext.mod.full)

bayesplot::pp_check(ext.mod.full)

plot(ext.mod.full$fitted.values)

posts = as.data.frame(ext.mod.full)

# NDD in large, high diversity
mean(posts$alphandd)

# NDD in large, low diversity
mean(posts$alphandd + posts$`low.varlow.div:alphandd`)

# NDD in small, high diversity populations
mean(posts$alphandd + posts$`n.pop0small:alphandd`)

ndd.posts = with(
  posts,
  cbind(
    large.hivar = alphandd,
    large.lovar = alphandd + `low.varlow.div:alphandd`,
    small.hivar = alphandd + `n.pop0small:alphandd`,
    small.lovar = alphandd + `low.varlow.div:alphandd` + `n.pop0small:alphandd` + `n.pop0small:low.varlow.div:alphandd`
    )
  )

apply(ndd.posts, 2, mean) %>% round(2)
apply(ndd.posts, 2, function(x) quantile(x, c(0.025, 0.975))) %>% round(2)

size.posts = with(
  posts,
  cbind(
    hivar.dindp = n.pop0small,
    hivar.ddepn = n.pop0small + `n.pop0small:alphandd`,
    lovar.dindp = n.pop0small + `n.pop0small:low.varlow.div`,
    lovar.ddepn = n.pop0small + `n.pop0small:alphandd` + `n.pop0small:low.varlow.div` + `n.pop0small:low.varlow.div:alphandd`
  )
)

apply(size.posts, 2, mean) %>% round(2)
apply(size.posts, 2, function(x) quantile(x, c(0.025, 0.975))) %>% round(2)

gdiv.posts = with(
  posts,
  cbind(
    large.dindp = low.varlow.div,
    large.ddepn = low.varlow.div + `low.varlow.div:alphandd`,
    small.dindp = low.varlow.div + `n.pop0small:low.varlow.div`,
    small.ddepn = low.varlow.div + `low.varlow.div:alphandd` + `n.pop0small:low.varlow.div` + `n.pop0small:low.varlow.div:alphandd`
  )
)

apply(gdiv.posts, 2, mean) %>% round(2)
apply(gdiv.posts, 2, function(x) quantile(x, c(0.025, 0.975))) %>% round(2)

### Get estimates for 

geno.ext.fulls = expand.grid(gbar = (-5:5)/10,
                             alpha = c("di", "ndd"),
                             n.pop0 = c("large", "small"),
                             low.var = c("low.div", "hi.div"))

epreds = posterior_epred(ext.mod.full, newdata = geno.ext.fulls,
                         seed = 196200, draws = 200) %>%
  as.data.frame() %>%
  t() %>%
  cbind(geno.ext.fulls, .)

epreds = epreds %>%
  gather(key = draw, value = estimate, -c(gbar, alpha, n.pop0, low.var))

geno.full.row = ggplot(epreds, aes(x = 2.75 - gbar, y = estimate)) +
  geom_line(
    aes(
      group = interaction(alpha, draw),
      colour = factor(alpha)
    ),
    size = 0.1
  ) +
  geom_rug(
    data = ext.data,
    inherit.aes = FALSE,
    aes(
      x = 2.75 - gbar
    ),
    sides = "b",
    size = 0.1
  ) +
  labs(x = 'Initial maladaptation', y = '') +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        plot.margin = margin(b = 0, r = 5, l = 20, unit = 'pt'))

geno.full.row

data.plots.row = plot_grid(size.plot.row, inst.plot.row, 
                           cuml.plot.row, geno.full.row, 
                           labels = c('(A)', '(B)', '(C)', '(D)'),
                           nrow = 4)

data.plots.row

plot_grid(data.plots.row, extinct.legend, ncol = 1, rel_heights = c(1, .1)) %>%
  save_plot(filename = 'simulations/analysis_results/figure_drafts/draft_figs/fig_ext_fig.png',
            base_width = 8, base_height = 8)

gbar.slopes = with(
  -1*posts,
  cbind(
    large.hivar.did = gbar,
    small.hivar.did = gbar + `n.pop0small:gbar`,
    large.lovar.did = gbar + `low.varlow.div:gbar`,
    large.hivar.ndd = gbar + `alphandd:gbar`,
    small.lovar.did = gbar + `n.pop0small:gbar` + `low.varlow.div:gbar` + `n.pop0small:low.varlow.div:gbar`,
    small.hivar.ndd = gbar + `n.pop0small:gbar` + `alphandd:gbar` + `n.pop0small:alphandd:gbar`,
    large.lovar.ndd = gbar + `low.varlow.div:gbar` + `alphandd:gbar` + `low.varlow.div:alphandd:gbar`,
    small.lovar.ndd = gbar + `n.pop0small:gbar` + `low.varlow.div:gbar` + `alphandd:gbar` +
                      `n.pop0small:low.varlow.div:gbar` + `n.pop0small:alphandd:gbar` + `low.varlow.div:alphandd:gbar` + 
                      `n.pop0small:low.varlow.div:alphandd:gbar`
  )
)

apply(gbar.slopes, 2, mean) %>% round(2)
apply(gbar.slopes, 2, function(x) quantile(x, c(0.025, 0.975))) %>% round(2)
