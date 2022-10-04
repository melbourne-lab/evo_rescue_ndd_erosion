# Script containing all results pertaining to extinctions in the long (50
# generation) simulations.
# Contains code for Figure 2
# SN - 29 Mar 2021
# re-run for archiving 4 Oct 2022

# Load packages
library(rstanarm)
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
    read.csv('simulations/outputs/longsims_n100_a000_hivar.csv'),
    read.csv('simulations/outputs/longsims_n100_a000_lowvar.csv'),
    read.csv('simulations/outputs/longsims_n20_a000_hivar.csv'),
    read.csv('simulations/outputs/longsims_n20_a000_lowvar.csv'),
    read.csv('simulations/outputs/longsims_n100_a035_hivar.csv'),
    read.csv('simulations/outputs/longsims_n100_a035_lowvar.csv'),
    read.csv('simulations/outputs/longsims_n20_a035_hivar.csv'),
    read.csv('simulations/outputs/longsims_n20_a035_lowvar.csv')
  ) %>%
  # Re-index time to start at generation 0 instead of 1
  mutate(gen = gen - 1)

# Add in extinction designation
all.data = all.data %>%
  # Get final generation's record for each population
  group_by(trial, n.pop0, low.var, alpha)  %>%
  slice_max(gen) %>%
  # Add "extinct" denotation
  group_by(n.pop0, low.var, alpha, trial) %>%
  summarise(extinct = gen < 50 & n < 5000) %>%
  # Combine this "extinct" column with main data frame
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
  # Reformating levels into plot-readible labels
  mutate(romnum = factor(ifelse(n.pop0 %in% 100,
                                ifelse(low.var, '(i)', '(ii)'),
                                ifelse(low.var, '(iii)', '(iv)'))),
         n.pop0 = factor(ifelse(n.pop0 %in% 100, "Large", "Small")),
         alpha = factor(ifelse(alpha > 0, "Density dependent", "Density independent")),
         low.var = factor(ifelse(low.var, "Low diversity", "High diversity")))

##### Bayesian model for Figure 3D and effect sizes

# Get data for extinction models:
# one row (entry) per simulation - using the first generation here to capture
# initial genotype.
ext.data = all.data %>% 
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  mutate(
    low.var = factor(ifelse(low.var, 'low.div', 'hi.div')),
    n.pop0  = factor(ifelse(n.pop0 > 50, 'large', 'small')),
    alpha   = factor(ifelse(alpha > 0, 'ndd', 'di'))
  )

# Fit model
ext.mod.full =  stan_glm(
  extinct ~ n.pop0 * low.var * alpha * gbar,
  data = ext.data, 
  seed = 992,
  cores = 4,
  prior = normal(0, 5),
  family = 'binomial'
)

summary(ext.mod.full)

# Model checks
bayesplot::pp_check(ext.mod.full)
plot(ext.mod.full$fitted.values)

# Posterior estimates for getting effect sizes
posts = as.data.frame(ext.mod.full)

### Effect sizes for Table 2

# Effect sizes for density dependence for all combos of other treatments

ndd.posts = with(
  posts,
  cbind(
    large.hivar = alphandd,
    large.lovar = alphandd + `low.varlow.div:alphandd`,
    small.hivar = alphandd + `n.pop0small:alphandd`,
    small.lovar = alphandd + `low.varlow.div:alphandd` + `n.pop0small:alphandd` + `n.pop0small:low.varlow.div:alphandd`
  )
)

# Mean and 95% CI
apply(ndd.posts, 2, mean) %>% round(2)
apply(ndd.posts, 2, function(x) quantile(x, c(0.025, 0.975))) %>% round(2)

# Effect sizes for size (small) for all combos of other treatments

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

# Effect sizes for low genetic diversity for all combos of other treatments

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

# Effect s izes for increasing maladaptation

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

# Arranged to highlight NDD effects
apply(gbar.slopes, 2, mean)[c(1,4,3,7,2,6,5,8)] %>% round(2)
apply(gbar.slopes, 2, function(x) quantile(x, c(0.025, 0.975)))[,c(1,4,3,7,2,6,5,8)] %>% round(2)

### Other summary statistics

# Increased extinction risk for each parameter combo
ext.data %>% 
  group_by(n.pop0, low.var, alpha) %>%
  summarise(pext = mean(extinct)) %>%
  spread(key = alpha, value = pext) %>%
  mutate(incr = ndd / di)

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
  geom_line(
    aes(
      y = p.extinct,
      group = factor(alpha),
      colour = factor(alpha)
    )
  ) +
  labs(x = 'Generation', y = '') +
  scale_fill_manual(values = c('purple', 'black')) +
  scale_colour_manual(values = c('purple', 'black')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  theme(legend.position = 'none',        
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        plot.margin = margin(b = 0, r = 5, l = 20, unit = 'pt'))

### Get estimates for Fig. 3D (posterior estimates of extinction probability
#   according to GLM)

# Get "genotypes" for which to get 
geno.ext.fulls = expand.grid(
  gbar = (-5:5)/10,
  alpha = c("di", "ndd"),
  n.pop0 = c("large", "small"),
  low.var = c("low.div", "hi.div")
) %>%
  filter(n.pop0 %in% 'small' | abs(gbar) < 0.25)

epreds = posterior_epred(
  ext.mod.full, 
  newdata = geno.ext.fulls,
  seed = 196200, 
  draws = 200
) %>%
  as.data.frame() %>%
  t() %>%
  cbind(geno.ext.fulls, .)

epreds = epreds %>%
  gather(key = draw, value = estimate, -c(gbar, alpha, n.pop0, low.var))

geno.full.row = ggplot(epreds, aes(x = 2.8 - gbar, y = estimate)) +
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
      x = 2.8 - gbar
    ),
    sides = "b",
    size = 0.1
  ) +
  labs(x = 'Initial maladaptation', y = '') +
  scale_color_manual(values = c('black', 'purple')) +
  scale_x_continuous(breaks = 2 * (12:16) / 10) +
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
                           label_x = c(-0.01, -0.01, -0.01, -0.01),
                           label_size = 12,
                           nrow = 4)

data.plots.row

# Get legend for figure
extinct.legend = get_legend(
  inst.plot.row + 
    guides(fill = guide_legend(''), colour = guide_legend('')) +
    theme(legend.position = 'bottom',
          legend.box.margin = margin(1, 0, 0, 0))
)

# Export plot
plot_grid(data.plots.row, extinct.legend, ncol = 1, rel_heights = c(1, .1)) %>%
  save_plot(filename = 'analysis_results/figures/fig_ext_fig.png',
            base_width = 8, base_height = 8)
