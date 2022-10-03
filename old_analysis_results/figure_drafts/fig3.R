### Genotypic change figure

# Load packages

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

# Load data for intrinsic fitness

ext.d = read.csv('simulations/outputs/final_results/alldata_combined.csv') %>%
  mutate(w = wbar) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct) %>%
  summarise(wbar = mean(w),
            wvar = var(w),
            n = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("High variance", "Low variance")))

# Plot fitness changes

ext.d %>%
  ggplot(aes(x = gen)) +
  geom_segment(
    aes(x = 0, xend = 15, y = 1, yend = 1),
    linetype = 3, colour = 'gray'
  ) +
  geom_line(
    aes(
      y = wbar,
      group = interaction(n.pop0, low.var, alpha, extinct),
      linetype = extinct,
      colour = alpha
    )
  ) +
  geom_ribbon(
    aes(
      ymin = wbar - 2 * sqrt(wvar / n),
      ymax = wbar + 2 * sqrt(wvar / n),
      group = interaction(n.pop0, low.var, alpha, extinct),
      fill = alpha
    ),
    alpha = 0.25
  ) +
  facet_wrap( ~ paste(reorder(n0, desc(n0)), low.var, sep = ', '), ncol = 2) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_linetype(labels = c('Surviving', 'Extinct')) +
  labs(x = 'Generation', y = expression(paste('Mean intrinsic fitness, ', W[t]))) +
  guides(colour = guide_legend('', nrow = 2), 
         fill = guide_legend('', nrow = 2), 
         linetype = guide_legend('', nrow = 2)) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        panel.grid = element_line(colour = 'gray88'),
        legend.position = 'bottom',
        strip.background = element_rect(colour = 'black')) +
  ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_3_w.pdf',
         height = 5, width = 5)

# Supplemental: get differences in extant variances

ext.d %>%
  filter(!extinct) %>%
  select(-c(wvar, n)) %>%
  spread(key = alpha, value = wbar) %>%
  rename(di = `Density independent`,
         dd = `Density dependent`) %>%
  mutate(wdiff = dd - di) %>%
  summarise(maxd = max(wdiff))

# Supplemental: first time period with wbar > 1
ext.d %>%
  filter(!extinct) %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(rescue.time = min(gen[wbar > 1]))

ext.g %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = gbar,
      group = interaction(n.pop0, low.var, alpha, extinct),
      linetype = extinct,
      colour = alpha
    )
  ) +
  geom_ribbon(
    aes(
      ymin = gbar - 2 * sqrt(gvar / n),
      ymax = gbar + 2 * sqrt(gvar / n),
      group = interaction(n.pop0, low.var, alpha, extinct),
      fill = alpha
    ),
    alpha = 0.25
  ) +
  facet_wrap( ~ paste(reorder(n0, desc(n0)), low.var, sep = ', '), ncol = 2) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_linetype(labels = c('Surviving', 'Extinct')) +
  labs(x = 'Generation', y = expression(paste('Mean intrinsic fitness, ', W[t]))) +
  guides(colour = guide_legend('', nrow = 2), 
         fill = guide_legend('', nrow = 2), 
         linetype = guide_legend('', nrow = 2)) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        panel.grid = element_line(colour = 'gray88'),
        legend.position = 'bottom',
        strip.background = element_rect(colour = 'black'))

### Same plots but now with genotype instead of intrinsic fitn.

ext.g = read.csv('simulations/outputs/final_results/alldata_combined.csv') %>%
  mutate(g = gbar) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct) %>%
  summarise(gbar = mean(g),
            gvar = var(g),
            n = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("High variance", "Low variance")))

ext.g %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = gbar,
      group = interaction(n.pop0, low.var, alpha, extinct),
      linetype = extinct,
      colour = alpha
    )
  ) +
  geom_ribbon(
    aes(
      ymin = gbar - 2 * sqrt(gvar / n),
      ymax = gbar + 2 * sqrt(gvar / n),
      group = interaction(n.pop0, low.var, alpha, extinct),
      fill = alpha
    ),
    alpha = 0.25
  ) +
  facet_wrap( ~ paste(reorder(n0, desc(n0)), low.var, sep = ', '), ncol = 2) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_linetype(labels = c('Surviving', 'Extinct')) +
  labs(x = 'Generation', y = expression(paste('Mean intrinsic fitness, ', W[t]))) +
  guides(colour = guide_legend('', nrow = 2), 
         fill = guide_legend('', nrow = 2), 
         linetype = guide_legend('', nrow = 2)) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        panel.grid = element_line(colour = 'gray88'),
        legend.position = 'bottom',
        strip.background = element_rect(colour = 'black')) +
  ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_3_g.pdf',
         height = 5, width = 5)

# Combined genotype over time and genotypic extinction plot

grid.g.res = 0.1

g.ext = all.data %>%
  # Now, bin variances and aggregate extinction probabilities
  mutate(gr = round((1/grid.g.res) * gbar) / (1/grid.g.res)) %>%
  group_by(gen, n.pop0, low.var, alpha, gr) %>%
  summarise(p.extinct = mean(extinct),
            n = n())

# Plot with no alpha (gives no sense of sample size of each grid cell)

extinct.plot1 = 
  g.ext %>%
  filter(gen < 15) %>%
  ggplot() +
  geom_tile(
    aes(
      x = gen,
      y = gr,
      fill = p.extinct
    )
  ) +
  facet_wrap(alpha ~ paste(n.pop0, low.var, sep = ', '), nrow = 2) +
  scale_fill_viridis_b(option = 'B') +
  labs(x = 'Generation', y = 'Mean population fitness') +
  guides(fill = guide_colorbar('Extinction\nprobability', barwidth = unit(3, 'cm'))) +
  theme(panel.background = element_rect(fill = NA),
        legend.position = 'bottom')

# Original plot

g.mean.plot = ext.g %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = gbar,
      group = interaction(n.pop0, low.var, alpha, extinct),
      linetype = extinct,
      colour = alpha
    )
  ) +
  geom_ribbon(
    aes(
      ymin = gbar - 2 * sqrt(gvar / n),
      ymax = gbar + 2 * sqrt(gvar / n),
      group = interaction(n.pop0, low.var, alpha, extinct),
      fill = alpha
    ),
    alpha = 0.25
  ) +
  facet_wrap( ~ paste(reorder(n0, desc(n0)), low.var, sep = ', '), ncol = 4) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_linetype(labels = c('Surviving', 'Extinct')) +
  labs(x = 'Generation', y = expression(paste('Mean intrinsic fitness, ', W[t]))) +
  guides(colour = guide_legend(''), 
         fill = guide_legend(''), 
         linetype = guide_legend('')) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        panel.grid = element_line(colour = 'gray88'),
        legend.position = 'bottom',
        strip.background = element_rect(colour = 'black'))

plot_grid(g.mean.plot, extinct.plot1, nrow = 2, rel_heights = c(1, 2))

