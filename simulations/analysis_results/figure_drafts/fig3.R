### Genotypic change figure

# Load packages

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

# Load data

ext.d = read.csv('simulations/outputs/final_results/alldata_combined.csv') %>%
  mutate(w = wbar) %>%
  group_by(n.pop0, low.var, alpha, gen, extinct) %>%
  summarise(wbar = mean(w),
            wvar = var(w),
            n = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("Low variance", "High variance")))

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
  scale_linetype(labels = c('Extant', 'Extinct')) +
  labs(x = 'Generation', y = expression(paste('Mean intrinsic fitness, ', W[t]))) +
  guides(colour = guide_legend('', nrow = 2), 
         fill = guide_legend('', nrow = 2), 
         linetype = guide_legend('', nrow = 2)) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        panel.grid = element_line(colour = 'gray88'),
        legend.position = 'bottom',
        strip.background = element_rect(colour = 'black')) +
  ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_3.pdf',
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
