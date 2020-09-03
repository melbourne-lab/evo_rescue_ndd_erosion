### Population size and dynamics figures.

# Load packages

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

# Load data

all.n = rbind(
  read.csv('simulations/outputs/final_results/a000_hivar_extinctions.csv') %>%
    mutate(low.var = FALSE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a000_lovar_extinctions.csv') %>%
    mutate(low.var = TRUE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a035_hivar_extinctions.csv') %>%
    mutate(low.var = FALSE, ndd = TRUE),
  read.csv('simulations/outputs/final_results/a035_lovar_extinctions.csv') %>%
    mutate(low.var = TRUE, ndd = TRUE)
) %>%
  group_by(low.var, ndd, n.pop0) %>%
  mutate(n.trials = pext * 8000,
         # Number of trials going extinct
         d.trials = -c(0, diff(n.trials)),
         # Proportion of extant trials going extinct
         p.ins1 = c(d.trials[-1], 0) / n.trials,
         # Proportion of all trials going extinct
         p.ins2 = c(d.trials[-1], 0) / 8000,
         # Proportion of extinct trials extinct in this generation
         p.ins3 = c(d.trials[-1], 0) / sum(d.trials)) %>%
 # filter(gen < 15) %>%
  ungroup() %>%
  mutate(init.size = factor(n.pop0, labels = c('Initially small', 'Initially large')),
         low.var = factor(low.var, labels = c('High variation', 'Low variation')),
         ndd = factor(ndd, labels = c('Density independent', 'Density dependent')))

all.n

# Try figure 1: population size over time

all.n %>%
  ggplot(aes(x = gen)) +
  geom_segment(
    aes(
      x = 0, xend = 15,
      y = n.pop0, yend = n.pop0
    ),
    linetype = 3,
    size = 0.5,
    colour = 'gray'
  ) +
  geom_line(
    aes(
      y = nbar,
      group = interaction(ndd, n.pop0, low.var),
      colour = factor(ndd),
      linetype = low.var
    )
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / n.trials),
      ymax = nbar + 2 * sqrt(nvar / n.trials),
      group = interaction(ndd, n.pop0, low.var),
      colour = factor(ndd),
      fill = factor(ndd)
    ),
    alpha = 0.2,
    size = 0.125
  ) +
  scale_color_manual(
    values = c('black', 'purple'),
    labels = c("Density\nindependent", "Density\ndependent"),
    name = "Growth form"
  ) +
  scale_fill_manual(
    values = c('black', 'purple'),
    labels = c("Density\nindependent", "Density\ndependent"),
    name = "Growth form"
  ) +
  scale_linetype(
    name = "Genetic variation",
    labels = c("High", "Low")
  ) +
  labs(x = 'Generation', y = 'Mean population size') +
  scale_y_log10() +
  facet_wrap(~ init.size) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = 'bottom') +
  ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_1.pdf',
         width = 8, height = 5)

### Extinction plot

inst.probs = ggplot(all.d %>% filter(gen < 15), aes(x = gen)) +
  geom_line(
    aes( 
      y = p.ins1, 
      group = interaction(n.pop0, bottleneck, ndd),
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = p.ins1 - 2 * sqrt(p.ins1 * (1-p.ins1) / n.trials),
      ymax = p.ins1 + 2 * sqrt(p.ins1 * (1-p.ins1) / n.trials),
      group = interaction(n.pop0, bottleneck, ndd),
      fill = ndd
    ),
    alpha = 0.2
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap( ~ paste(n.pop0, bottleneck, sep = ', '), ncol = 1) +
  labs(x = 'Generation', y = 'Probability of extinction') +
  theme(legend.position = 'none',
        panel.background = element_blank())

cuml.probs = ggplot(all.d %>% filter(gen < 15), aes(x = gen)) +
  geom_ribbon(
    aes(
      xmin = 0, xmax = 14,
      ymin = 0, ymax = 1 - pext,
      group = interaction(n.pop0, bottleneck, ndd),
      fill = ndd
    ),
    alpha = 0.25
  ) +
  scale_fill_manual(values = c('black', 'purple')) +
  labs(x = 'Generation') +
  facet_wrap( ~ paste(n.pop0, bottleneck, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.background = element_blank())

extinct.legend = get_legend(
  inst.probs + 
    guides(fill = guide_legend(''), colour = guide_legend('')) +
    theme(legend.position = 'bottom',
          legend.box.margin = margin(1, 0, 0, 0))
)

data.plots = plot_grid(inst.probs, cuml.probs)
plot_grid(data.plots, extinct.legend, ncol = 1, rel_heights = c(1, .1))

plot_grid(data.plots, extinct.legend, ncol = 1, rel_heights = c(1, .1)) %>%
  save_plot(filename = 'simulations/analysis_results/figure_drafts/draft_figs/fig_2.pdf',
            base_width = 5, base_height = 6)
