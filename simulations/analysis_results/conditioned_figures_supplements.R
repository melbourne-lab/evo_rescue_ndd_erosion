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

