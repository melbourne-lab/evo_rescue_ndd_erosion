### Plots associated with genetic diversity and fixation
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

### Read in genetic data

# Complete pooling data (all in one curve/average)
all.g = rbind(
  read.csv('simulations/outputs/final_results/a000_hivar_gene_all.csv')  %>%
    mutate(low.var = FALSE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a035_hivar_gene_all.csv') %>%
    mutate(low.var = FALSE, ndd = TRUE),
  read.csv('simulations/outputs/final_results/a000_lowvar_gene_all.csv') %>%
    mutate(low.var = TRUE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a035_lowvar_gene_all.csv') %>%
    mutate(low.var = TRUE, ndd = TRUE)
) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c('Initially small', 'Initially large')),
         low.var = factor(low.var, labels = c('High variation', 'Low variation')),
         ndd = factor(ndd, labels = c('Density independent', 'Density dependent')))

# Pooling by extinct/extant

ext.g = rbind(
  read.csv('simulations/outputs/final_results/a000_hivar_gene_ext.csv')  %>%
    mutate(low.var = FALSE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a035_hivar_gene_ext.csv') %>%
    mutate(low.var = FALSE, ndd = TRUE),
  read.csv('simulations/outputs/final_results/a000_lowvar_gene_ext.csv') %>%
    mutate(low.var = TRUE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a035_lowvar_gene_ext.csv') %>%
    mutate(low.var = TRUE, ndd = TRUE)
) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c('Initially small', 'Initially large')),
         low.var = factor(low.var, labels = c('High variation', 'Low variation')),
         ndd = factor(ndd, labels = c('Density independent', 'Density dependent')))

# Pooling by extant/extinct generation

gen.g = rbind(
  read.csv('simulations/outputs/final_results/a000_hivar_gene_gen.csv')  %>%
    mutate(low.var = FALSE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a035_hivar_gene_gen.csv') %>%
    mutate(low.var = FALSE, ndd = TRUE),
  read.csv('simulations/outputs/final_results/a000_lowvar_gene_gen.csv') %>%
    mutate(low.var = TRUE, ndd = FALSE),
  read.csv('simulations/outputs/final_results/a035_lowvar_gene_gen.csv') %>%
    mutate(low.var = TRUE, ndd = TRUE)
) %>%
  mutate(n.pop0 = factor(n.pop0, labels = c('Initially small', 'Initially large')),
         low.var = factor(low.var, labels = c('High variation', 'Low variation')),
         ndd = factor(ndd, labels = c('Density independent', 'Density dependent')))


### Genetic variation plot

comp.pool.plot = all.g %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = vbar,
      group = interaction(n.pop0, low.var, ndd),
      colour = ndd
    )
  ) +
  geom_ribbon(
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, low.var, ndd),
      fill = ndd
    ),
    alpha = 0.1
  ) +
  scale_y_continuous(limits = c(0, 0.5)) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_colour_manual(values = c('black', 'purple')) +
  guides(colour = guide_legend("Growth form", nrow = 1),
         fill = guide_legend("Growth form", nrow = 1)) +
  labs(x = '', y = 'Genetic variation') +
  facet_wrap( ~ paste(low.var, n.pop0, sep = ', '), ncol = 4) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(colour = 'gray88'),
        panel.border = element_rect(fill = NA),
        legend.position = 'none')

comp.pool.plot

gent.pool.plot = gen.g %>%
  ggplot(aes(x = gen)) +
  geom_line(
    data = . %>% filter(ext.gen < 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, low.var, ndd, ext.gen),
      colour = ext.gen
    )
  ) +
  geom_ribbon(
    data = . %>% filter(ext.gen < 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, low.var, ndd, ext.gen),
      fill = ext.gen
    ),
    alpha = 0.1
  ) +
  geom_line(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, low.var, ndd)
    ),
    colour = 'black'
  )  +
  geom_ribbon(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, low.var, ndd)
    ),
    fill = 'black',
    alpha = 0.1
  ) +
  scale_colour_gradient(low = 'lightpink1', high = 'darkred') +
  scale_fill_gradient(low = 'lightpink1', high = 'darkred') +
  guides(colour = guide_legend("Extinction\ngeneration"),
         fill = guide_legend('Extinction\ngeneration')) +
  labs(x = 'Generation', y = 'Genetic variation') +
  facet_wrap(ndd ~ paste(low.var, n.pop0, sep = ', '), ncol = 4) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid = element_line(colour = 'gray88'),
        legend.position = 'none')

gent.pool.plot

comb.plots = plot_grid(comp.pool.plot, NULL, gent.pool.plot, 
                       ncol = 1, rel_heights = c(1, -0.1, 2),
                       labels = c('(A)', '', '(B)'))

comb.plots %>%
  save_plot(filename = 'simulations/analysis_results/figure_drafts/draft_figs/fig_4.pdf',
            base_height = 8, base_width = 8)

### Fixation plot

# Fixation of only extant populations

ext.g %>%
  filter(!extinct) %>%
  gather(key = fixtype, value = p, p.fix.pos, p.fix.neg) %>%
  mutate(fixtype = factor(fixtype, labels = c('Negative allele', 'Positive allele'))) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = p,
      group = interaction(n.pop0, low.var, ndd, fixtype),
      linetype = paste(n.pop0, low.var, sep = ', '),
      colour = ndd
    )
  ) +
  scale_linetype_manual(values = c(4, 6, 2, 5)) +
  scale_colour_manual(values = c('black', 'purple')) +
  facet_wrap(~ reorder(fixtype, desc(fixtype))) +
  guides(linetype = guide_legend('', nrow = 2),
         colour = guide_legend('', nrow = 2)) +
  labs(x = 'Generation', y = 'Probability of fixation') +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.direction = 'horizontal',
        legend.position = 'bottom',
        legend.text = element_text(size = 8),
        strip.background = element_rect(colour = 'black'),
        strip.text = element_text(size = 12)) +
  ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_5_extant.pdf',
         height = 5, width = 8)

all.g %>%
  gather(key = fixtype, value = p, p.fix.pos, p.fix.neg) %>%
  mutate(fixtype = factor(fixtype, labels = c('Negative allele', 'Positive allele'))) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = p,
      group = interaction(n.pop0, low.var, ndd, fixtype),
      linetype = paste(n.pop0, low.var, sep = ', '),
      colour = ndd
    )
  ) +
  # This isn't a binomial random variable... or is it?? these intervals look off.
  # geom_ribbon(
  #   aes(
  #     ymin = p - 2 * sqrt(p * (1 - p) / n),
  #     ymax = p + 2 * sqrt(p * (1 - p) / n),
  #     group = interaction(n.pop0, low.var, ndd, fixtype),
  #     linetype = paste(n.pop0, low.var, sep = ', '),
  #     fill = ndd
  #   ),
  #   alpha = 0.2
  # ) +
  scale_linetype_manual(values = c(4, 6, 2, 5)) +
  scale_colour_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap(~ reorder(fixtype, desc(fixtype))) +
  guides(linetype = guide_legend('', nrow = 2),
         colour = guide_legend('', nrow = 2),
         fill = guide_legend('', nrow = 2)) +
  labs(x = 'Generation', y = 'Probability of fixation') +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.direction = 'horizontal',
        legend.position = 'bottom',
        legend.text = element_text(size = 8),
        strip.background = element_rect(colour = 'black'),
        strip.text = element_text(size = 12)) +
  ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_5.pdf',
         height = 5, width = 8)

