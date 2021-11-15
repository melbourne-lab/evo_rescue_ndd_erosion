# Supplemental figures for various figures of:
#   - state variables conditioned on extinction/survival at end of simulation
#   - state variables conditioned on survival/generation of extinction

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

# Load in all simulation data
all.data = read.csv('simulations/outputs/final_results/alldata_combined.csv')

### Fig C1: population size

# Summarise data by param combo, generation, extinction status
# n.b. need to add zeros for extinct populations
ext.n = all.data %>%
  # Get relevant columns only
  select(trial, gen, n, n.pop0, low.var, alpha) %>%
  # Add extinctions in (for mean population size)
  rbind(
    expand.grid(trial = 1:4000,
                gen = 1:15,
                n = 0,
                n.pop0 = c(20, 100),
                low.var = c(TRUE, FALSE),
                alpha = c(0, 0.0035))
  ) %>%
  group_by(trial, gen, n.pop0, low.var, alpha) %>%
  summarise(n = sum(n)) %>%
  # Get rid of trials where population was erroneously added
  group_by(trial, n.pop0, low.var, alpha) %>%
  filter(any(n>0)) %>%
  # Add extinction status to trial
  group_by(trial, n.pop0, low.var, alpha) %>%
  mutate(ext = any(!n)) %>%
  ungroup() %>%
  # Aggregate (mean and variance in each generation)
  group_by(gen, n.pop0, low.var, alpha, ext) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            n.trials = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Initially small", "Initially large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("Low genetic diversity", "High genetic diversity"))) %>%
  # For purposes of plotting on log scale, remove records where nbar is 0
  filter(nbar > 0)

ext.n %>%
  ggplot(aes(x = gen, y = nbar)) +
  geom_line(
    aes(
      colour = factor(alpha),
      group = interaction(ext, alpha),
      linetype = ext
      )
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000),
                labels = c(0.1, 1, 10, 100, 1000)) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.grid.minor = element_blank(),
        legend.position = 'none') #,
        # strip.text = element_blank(),
        # strip.background = element_blank())

### Conditioning mean genotype on extinction

ext.gz = all.data %>%
  rename(g = gbar,
         z = zbar) %>%
  select(trial, gen, g, z, n.pop0, alpha, low.var, extinct) %>%
  gather(key = gz, value = val, c(g, z)) %>%
  # Aggregate (mean and variance in each generation)
  group_by(gen, n.pop0, low.var, alpha, extinct, gz) %>%
  summarise(varbar = mean(val),
            varvar = var(val),
            n.trials = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("High diversity", "Low diversity")))

ext.gz %>%
  ggplot(aes(x = gen, y = varbar)) +
  geom_ribbon(
    aes(
      ymin = varbar - 2 * sqrt(varvar / n.trials),
      ymax = varbar + 2 * sqrt(varvar / n.trials),
      fill = factor(alpha),
      group = interaction(extinct, alpha)
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      colour = factor(alpha),
      group = interaction(extinct, alpha),
      linetype = extinct
    )
  ) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap(gz ~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.grid.minor = element_blank(),
        legend.position = 'none')

# A cowplot would probably work better than this

ext.g.plot = ext.gz %>%
  filter(gz %in% 'g') %>%
  ggplot(aes(x = gen, y = varbar)) +
  geom_ribbon(
    aes(
      ymin = varbar - 2 * sqrt(varvar / n.trials),
      ymax = varbar + 2 * sqrt(varvar / n.trials),
      fill = factor(alpha),
      group = interaction(extinct, alpha)
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      colour = factor(alpha),
      group = interaction(extinct, alpha),
      linetype = extinct
    ),
    size = 1.5
  ) +
  labs(x = '', y = expression('Mean population genotype, ' ~ bar(g[t]))) +
  scale_linetype_manual(values = c(1, 2)) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap( ~ paste(n0, low.var, sep = ', '), ncol = 4) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# ext.g.plot

ext.z.plot = ext.gz %>%
  filter(gz %in% 'z') %>%
  mutate(extinct = factor(extinct, labels = c('Surviving', 'Extinct'))) %>%
  ggplot(aes(x = gen, y = varbar)) +
  geom_ribbon(
    aes(
      ymin = varbar - 2 * sqrt(varvar / n.trials),
      ymax = varbar + 2 * sqrt(varvar / n.trials),
      fill = factor(alpha),
      group = interaction(extinct, alpha)
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      colour = factor(alpha),
      group = interaction(extinct, alpha),
      linetype = extinct
    ),
    size = 1.5
  ) +
  labs(x = 'Generation', y = expression('Mean population phenotype, ' ~ bar(z[t]))) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_linetype_manual(values = c(1, 2)) +
  facet_wrap( ~ paste(n0, low.var, sep = ', '), ncol = 4) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        strip.text = element_blank(),
        strip.background = element_blank())

# ext.z.plot

ext.gz.legend = get_legend(
  ext.z.plot + 
    guides(fill = guide_legend(''), colour = guide_legend(''),
           linetype = guide_legend(' ', values = c(1, 3))) +
    theme(legend.position = 'bottom',
          legend.box.margin = margin(1, 0, 0, 0))
)
# note: see and run "varn.legend" in fig4_fig5.R
# as this plot won't let us change the order or style of legends
# (truly a mystery)

ext.gz.plots = plot_grid(ext.g.plot, ext.z.plot, 
                        labels = c('(A)', '(B)'),
                        label_y = c(1, 1.075),
                        nrow = 2)
save_plot(
  plot_grid(ext.gz.plots, varn.legend, 
            ncol = 1, rel_heights = c(1, .1)),
  filename = 'simulations/analysis_results/figure_drafts/draft_figs/fig_supp_ext_gz.pdf',
  base_width = 8, base_height = 6
)

### Conditioning mean intrinsic fitness W on extinction

ext.wr = all.data %>%
  rename(w = wbar) %>%
  mutate(r = w * exp(-n * alpha)) %>%
  select(trial, gen, n.pop0, low.var, alpha, extinct, w, r) %>%
  gather(key = wr, value = val, c(w, r)) %>%
  # Aggregate (mean and variance in each generation)
  group_by(gen, n.pop0, low.var, alpha, extinct, wr) %>%
  summarise(varbar = mean(val),
            varvar = var(val),
            n.trials = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("High diversity", "Low diversity")))

ext.w.plot = ext.wr %>%
  filter(wr %in% 'w') %>%
  ggplot(aes(x = gen, y = varbar)) +
  geom_ribbon(
    aes(
      ymin = varbar - 2 * sqrt(varvar / n.trials),
      ymax = varbar + 2 * sqrt(varvar / n.trials),
      fill = factor(alpha),
      group = interaction(extinct, alpha)
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      colour = factor(alpha),
      group = interaction(extinct, alpha),
      linetype = extinct
    ),
    size = 1.5
  ) +
  scale_y_continuous(breaks = c(0.8, 1.0, 1.2, 1.4)) +
  labs(x = '', y = expression(atop('Mean population intrinsic', 'fitness, ' ~ bar(W[t])))) +
  scale_linetype_manual(values = c(1, 2)) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap( ~ paste(n0, low.var, sep = ', '), ncol = 4) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ext.r.plot = ext.wr %>%
  filter(wr %in% 'r') %>%
  mutate(extinct = factor(extinct, labels = c('Surviving', 'Extinct'))) %>%
  ggplot(aes(x = gen, y = varbar)) +
  geom_ribbon(
    aes(
      ymin = varbar - 2 * sqrt(varvar / n.trials),
      ymax = varbar + 2 * sqrt(varvar / n.trials),
      fill = factor(alpha),
      group = interaction(extinct, alpha)
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      colour = factor(alpha),
      group = interaction(extinct, alpha),
      linetype = extinct
    ),
    size = 1.5
  ) +
  labs(x = 'Generation', y = expression(atop('Mean population extrinsic', 'fitness, '~ bar(R[t])))) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_linetype_manual(values = c(1, 2)) +
  facet_wrap( ~ paste(n0, low.var, sep = ', '), ncol = 4) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        strip.text = element_blank(),
        strip.background = element_blank())

ext.wr.plots = plot_grid(ext.w.plot, ext.r.plot, 
                         labels = c('(A)', '(B)'),
                         label_y = c(1, 1.05),
                         nrow = 2)
save_plot(
  plot_grid(ext.wr.plots, varn.legend, 
            ncol = 1, rel_heights = c(1, .1)),
  filename = 'simulations/analysis_results/figure_drafts/draft_figs/fig_supp_ext_wr.pdf',
  base_width = 8, base_height = 6
)

### Conditioning on survival/*generation* of extinction

# Correct generation of extinction
# Define the generation of extinction to be the one where population size is 0/1
# If a population's final census is greater than 1, add another generation
all.data = all.data %>%
  group_by(trial, n.pop0, low.var, alpha) %>%
  mutate(ext.gen = max(gen) + as.numeric(n[which.max(gen)] > 1)) %>%
  mutate(extinct = !(max(gen) > 14 & all(n > 1)))

all.data %>%
  arrange(trial, n.pop0, low.var, alpha, desc(gen)) %>%
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE)

### Aggregate

gen.cond = all.data %>%
  group_by(n.pop0, low.var, alpha, ext.gen, extinct, gen) %>%
  rename(g = gbar, w = wbar) %>%
  summarise(
    # Population size
    nbar = mean(n),
    nvar = var(n),
    # Genotype
    gbar = mean(g),
    gvar = var(g),
    # Intrinsic fitness
    wbar = mean(w),
    wvar = var(w),
    # Fixation probabilities
    fixp = mean(p.fix.pos), # positive
    fixn = mean(p.fix.neg), # negative
    # Genetic variation
    vbar = mean(v),
    vvar = var(v),
    # Number of observations
    nobs = n()
  ) %>%
  ungroup()

head(gen.cond)

### Population size plot

gen.cond %>%
  mutate(
    n0 = factor(n.pop0, labels = c('Small', 'Large')),
    lv = factor(low.var, labels = c('High diversity', 'Low diversity')),
    dd = factor(alpha, labels = c('Density independent', 'Density dependent')),
    extinct = factor(extinct, labels = c("Surviving", "Extinct"))
  ) %>%
  ggplot(aes(x = gen)) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / nobs),
      ymax = nbar + 2 * sqrt(nvar / nobs),
      group = interaction(ext.gen, extinct),
      fill = extinct
    ),
    alpha = 0.1
  ) +
  geom_line(
    aes(
      y = nbar,
      group = interaction(ext.gen, extinct),
      colour = extinct
    )
  ) +
  scale_fill_manual(values = c('black', 'red'), '') +
  scale_color_manual(values = c('black', 'red'), '') +
  scale_y_log10() +
  labs(x = "Generation", y = "Population size") +
  facet_wrap(dd ~ paste(n0, lv, sep = ', '), ncol = 4) +
  theme(
    legend.position = 'bottom',
    panel.grid.major = element_line(colour = 'gray88'),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour = 'black')
  )

### Genotype

gen.cond %>%
  mutate(
    n0 = factor(n.pop0, labels = c('Small', 'Large')),
    lv = factor(low.var, labels = c('High diversity', 'Low diversity')),
    dd = factor(alpha, labels = c('Density independent', 'Density dependent')),
    extinct = factor(extinct, labels = c("Surviving", "Extinct"))
  ) %>%
  ggplot(aes(x = gen)) +
  # geom_ribbon(
  #   aes(
  #     ymin = gbar - 2 * sqrt(gvar / nobs),
  #     ymax = gbar + 2 * sqrt(gvar / nobs),
  #     group = interaction(ext.gen, extinct),
  #     fill = extinct
  #   ),
  #   alpha = 0.1
  # ) +
  geom_line(
    aes(
      y = gbar,
      group = interaction(ext.gen, extinct),
      colour = extinct
    )
  ) +
  scale_fill_manual(values = c('black', 'red'), '') +
  scale_color_manual(values = c('black', 'red'), '') +
  labs(x = "Generation", y = "Genotypee") +
  facet_wrap(dd ~ paste(n0, lv, sep = ', '), ncol = 4) +
  theme(
    legend.position = 'bottom',
    panel.grid.major = element_line(colour = 'gray88'),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour = 'black')
  )

### Intrinsic fitness (W)

gen.cond %>%
  mutate(
    n0 = factor(n.pop0, labels = c('Small', 'Large')),
    lv = factor(low.var, labels = c('High diversity', 'Low diversity')),
    dd = factor(alpha, labels = c('Density independent', 'Density dependent')),
    extinct = factor(extinct, labels = c("Surviving", "Extinct"))
  ) %>%
  ggplot(aes(x = gen)) +
  # geom_ribbon(
  #   aes(
  #     ymin = wbar - 2 * sqrt(wvar / nobs),
  #     ymax = wbar + 2 * sqrt(wvar / nobs),
  #     group = interaction(ext.gen, extinct),
  #     fill = extinct
  #   ),
  #   alpha = 0.1
  # ) +
  geom_line(
    aes(
      y = wbar,
      group = interaction(ext.gen, extinct),
      colour = extinct
    )
  ) +
  scale_fill_manual(values = c('black', 'red'), '') +
  scale_color_manual(values = c('black', 'red'), '') +
  labs(x = "Generation", y = "Intrinsic fitness") +
  facet_wrap(dd ~ paste(n0, lv, sep = ', '), ncol = 4) +
  theme(
    legend.position = 'bottom',
    panel.grid.major = element_line(colour = 'gray88'),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour = 'black')
  )

### Genetic variation

gen.cond %>%
  mutate(
    n0 = factor(n.pop0, labels = c('Small', 'Large')),
    lv = factor(low.var, labels = c('High diversity', 'Low diversity')),
    dd = factor(alpha, labels = c('Density independent', 'Density dependent')),
    extinct = factor(extinct, labels = c("Surviving", "Extinct"))
  ) %>% 
  ggplot(aes(x = gen)) +
  geom_ribbon(
    aes(
      ymin = vbar - 2 * sqrt(vvar / nobs),
      ymax = vbar + 2 * sqrt(vvar / nobs),
      group = interaction(ext.gen, extinct),
      fill = extinct
    ),
    alpha = 0.1
  ) +
  geom_line(
    aes(
      y = vbar,
      group = interaction(ext.gen, extinct),
      colour = extinct
    )
  ) +
  scale_fill_manual(values = c('black', 'red'), '') +
  scale_color_manual(values = c('black', 'red'), '') +
  labs(x = "Generation", y = "Genetic variation") +
  facet_wrap(dd ~ paste(n0, lv, sep = ', '), ncol = 4) +
  theme(
    legend.position = 'bottom',
    panel.grid.major = element_line(colour = 'gray88'),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour = 'black')
  )

### Fixation of alleles (positive)

gen.cond %>%
  mutate(
    n0 = factor(n.pop0, labels = c('Small', 'Large')),
    lv = factor(low.var, labels = c('High diversity', 'Low diversity')),
    dd = factor(alpha, labels = c('Density independent', 'Density dependent')),
    extinct = factor(extinct, labels = c("Surviving", "Extinct"))
  ) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = fixp,
      group = interaction(ext.gen, extinct),
      colour = extinct
    )
  ) +
  scale_fill_manual(values = c('black', 'red'), '') +
  scale_color_manual(values = c('black', 'red'), '') +
  labs(x = "Generation", y = "Loci at fixation (positive allele)") +
  facet_wrap(dd ~ paste(n0, lv, sep = ', '), ncol = 4) +
  theme(
    legend.position = 'bottom',
    panel.grid.major = element_line(colour = 'gray88'),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour = 'black')
  )

### Fixation of alleles (negative)

gen.cond %>%
  mutate(
    n0 = factor(n.pop0, labels = c('Small', 'Large')),
    lv = factor(low.var, labels = c('High diversity', 'Low diversity')),
    dd = factor(alpha, labels = c('Density independent', 'Density dependent')),
    extinct = factor(extinct, labels = c("Surviving", "Extinct"))
  ) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = fixn,
      group = interaction(ext.gen, extinct),
      colour = extinct
    )
  ) +
  scale_fill_manual(values = c('black', 'red'), '') +
  scale_color_manual(values = c('black', 'red'), '') +
  labs(x = "Generation", y = "Loci at fixation (negative allele)") +
  facet_wrap(dd ~ paste(n0, lv, sep = ', '), ncol = 4) +
  theme(
    legend.position = 'bottom',
    panel.grid.major = element_line(colour = 'gray88'),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour = 'black')
  )
