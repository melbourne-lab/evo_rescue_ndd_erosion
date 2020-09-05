### Genotypic change figure

# Load packages

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

# Load data

ext.d = rbind(
  read.csv('simulations/outputs/final_results/a000_hivar_type_ext.csv') %>%
    mutate(low.var = FALSE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a000_lovar_type_ext.csv') %>%
    mutate(low.var = TRUE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a035_hivar_type_ext.csv') %>%
    mutate(low.var = FALSE, ndd = TRUE),
  read.csv('simulations/outputs/final_results/a035_lovar_type_ext.csv') %>%
    mutate(low.var = TRUE, ndd = TRUE)
) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c('Initially small', 'Initially large')),
         low.var = factor(low.var, labels = c('High variation', 'Low variation')),
         ndd = factor(ndd, labels = c('Density independent', 'Density dependent')))

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
      group = interaction(n.pop0, low.var, ndd, extinct),
      linetype = extinct,
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = wbar - 2 * sqrt(wvar / n),
      ymax = wbar + 2 * sqrt(wvar / n),
      group = interaction(n.pop0, low.var, ndd, extinct),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  facet_wrap( ~ paste(low.var, reorder(n.pop0, desc(n.pop0)), sep = ', '), ncol = 2) +
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
