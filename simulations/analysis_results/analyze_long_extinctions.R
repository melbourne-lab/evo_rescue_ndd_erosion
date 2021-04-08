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
  mutate(n.pop0 = factor(ifelse(n.pop0 %in% 100, "Large", "Small")),
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

summary(ext.mod.final)

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
    size = 0.5) +
  labs(x = 'Generation', y = '') +
  scale_y_log10() +
  scale_alpha_manual(values = c(0, 1)) +
  scale_fill_manual(values = c('purple', 'black')) +
  scale_color_manual(values = c('purple', 'black')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  theme(legend.position = 'none',
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        strip.background = element_rect(colour = 'black'))

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
        panel.background = element_rect(fill = 'white'))

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
  labs(x = 'Generation', y = 'Probability of extinction') +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  theme(legend.position = 'none',        
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'))

geno.plot.row = geno.ext.preds %>%
  ggplot() +
  geom_line(
    aes(
      x = gbar, y = p, group = factor(alpha),
      colour = factor(alpha)
    )
  ) +
  labs(x = 'Initial genotype', y = '') +
  scale_color_manual(values = c('purple', 'black')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'))

data.plots.row = plot_grid(size.plot.row, cuml.plot.row, 
                           inst.plot.row, geno.plot.row, 
                           labels = c('(A)', '(B)', '(C)', '(D)'),
                           nrow = 4)

plot_grid(data.plots.row, extinct.legend, ncol = 1, rel_heights = c(1, .1)) %>%
  save_plot(filename = 'simulations/analysis_results/figure_drafts/draft_figs/fig_2_long_rows.pdf',
            base_width = 8, base_height = 8)
