### Plots associated with genetic diversity and fixation
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

### Read in genetic data

# All data
all.data = read.csv('simulations/outputs/final_results/alldata_combined.csv')

# Complete pooling data (all in one curve/average)
all.v = all.data %>%
  group_by(n.pop0, low.var, alpha, gen) %>%
  summarise(vbar = mean(v),
            vvar = var(v),
            n = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density dependent", "Density independent")),
         low.var = factor(low.var, labels = c("High diversity", "Low diversity")))

# Pooling by extinct/extant

ext.v = all.data %>%
  group_by(n.pop0, low.var, alpha, gen, extinct) %>%
  summarise(vbar = mean(v),
            vvar = var(v),
            n = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("High diversity", "Low diversity")))


# Pooling by extant/extinct generation

gen.v = all.data %>%
  group_by(n.pop0, low.var, alpha, gen, ext.gen) %>%
  summarise(vbar = mean(v),
            vvar = var(v),
            n = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("High diversity", "Low diversity")))

# Fixation probabilities, unconditional

all.fix = all.data %>%
  group_by(n.pop0, low.var, alpha, gen) %>%
  summarise(p.fix.pos = mean(p.fix.pos),
            p.fix.neg = mean(p.fix.neg),
            n = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("High diversity", "Low diversity")))

# Fixation probabilities by extinction

ext.fix = all.data %>%
  group_by(n.pop0, low.var, alpha, gen, extinct) %>%
  summarise(p.pos = mean(p.fix.pos),
            v.pos = var(p.fix.pos),
            p.neg = mean(p.fix.neg),
            v.neg = var(p.fix.neg),
            n = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("High diversity", "Low diversity")))

### Genetic variation plot

comp.pool.plot = all.v %>%
  mutate(alpha = relevel(alpha, "Density independent")) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = vbar,
      group = interaction(n.pop0, low.var, alpha),
      colour = alpha
    ),
    size = 1.25
  ) +
  geom_ribbon(
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, low.var, alpha),
      fill = alpha
    ),
    alpha = 0.1
  ) +
  scale_y_continuous(limits = c(0, 0.5)) +
  scale_fill_manual(values = c('purple', 'black')) +
  scale_colour_manual(values = c('purple', 'black')) +
  guides(colour = guide_legend("Growth form", nrow = 1),
         fill = guide_legend("Growth form", nrow = 1)) +
  labs(x = '', y = 'Genetic variance') +
  facet_wrap( ~ paste(n0, low.var, sep = ', '), ncol = 4) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(colour = 'gray88'),
        panel.border = element_rect(fill = NA),
        legend.position = 'none')

comp.pool.plot

gent.pool.plot = gen.v %>%
  ggplot(aes(x = gen)) +
  geom_line(
    data = . %>% filter(ext.gen < 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, low.var, alpha, ext.gen)
    ),
    colour = 'firebrick3'
  ) +
  geom_ribbon(
    data = . %>% filter(ext.gen < 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, low.var, alpha, ext.gen)
    ),
    fill = 'firebrick3',
    alpha = 0.1
  ) +
  geom_line(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      y = vbar,
      group = interaction(n.pop0, low.var, alpha)
    ),
    colour = 'black'
  )  +
  geom_ribbon(
    data = . %>% filter(ext.gen %in% 15),
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, low.var, alpha)
    ),
    fill = 'black',
    alpha = 0.1
  ) +
  # scale_colour_gradient(low = 'lightpink1', high = 'darkred') +
  # scale_fill_gradient(low = 'lightpink1', high = 'darkred') +
  guides(colour = guide_legend("Extinction\ngeneration"),
         fill = guide_legend('Extinction\ngeneration')) +
  labs(x = 'Generation', y = 'Genetic variance') +
  facet_wrap(alpha ~ paste(n0, low.var, sep = ', '), ncol = 4) +
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

# Figure 5

ext.fix %>%
  gather(key = fixtype.mean, value = p, p.pos, p.neg) %>%
  gather(key = fixtype.vars, value = v, v.pos, v.neg) %>%
  filter((grepl('pos', fixtype.mean) & grepl('pos', fixtype.vars)) |
           (grepl('neg', fixtype.mean) & grepl('neg', fixtype.vars)) ) %>%
  mutate(fixtype = gsub('^p\\.', '', fixtype.mean)) %>%
  select(-c(fixtype.mean, fixtype.vars)) %>%
  mutate(fixtype = factor(fixtype, labels = c('Maladaptive allele', 'Adaptive allele')),
         alpha = relevel(alpha, "Density dependent")) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = p,
      group = interaction(n.pop0, low.var, alpha, fixtype, extinct),
      linetype = extinct,
      colour = alpha
    ),
    size = 1.5
  ) +
  geom_ribbon(
    aes(
      ymin = p - 2 * sqrt(v / n),
      ymax = p + 2 * sqrt(v / n),
      group = interaction(n.pop0, low.var, alpha, fixtype, extinct),
      fill = alpha
    ),
    alpha = 0.2
  ) +
  scale_colour_manual(values = c('purple', 'black')) +
  scale_fill_manual(values = c('purple', 'black')) +
  scale_linetype(labels = c("Surviving", "Extinct")) +
  facet_wrap(reorder(fixtype, desc(fixtype)) ~ paste(n0, low.var, sep = ', '), ncol = 4) +
  guides(linetype = guide_legend('', nrow = 1),
         colour = guide_legend('', nrow = 1),
         fill = guide_legend('', nrow = 1)) +
  labs(x = 'Generation', y = 'Proportion of loci at fixation') +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.direction = 'horizontal',
        legend.position = 'bottom',
        legend.text = element_text(size = 12),
        strip.background = element_rect(colour = 'black'),
        strip.text = element_text(size = 12)) +
  ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_5_eight_line.pdf',
         width = 8, height = 5)

### Combined figure 4 where panel (a) is genetic variation, (b) is fixation rates

variance.plot = ext.v %>%
  mutate(alpha = relevel(alpha, "Density independent"),
         extinct = ifelse(extinct, 'Extinct', 'Surviving')) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = vbar,
      group = interaction(n.pop0, low.var, alpha, extinct),
      colour = alpha,
      linetype = extinct
    ),
    size = 1.25
  ) +
  geom_ribbon(
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, low.var, alpha, extinct),
      fill = alpha
    ),
    alpha = 0.1
  ) +
  scale_y_continuous(limits = c(0, 0.5)) +
  scale_colour_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_linetype_manual(values = c(2, 1)) +
  labs(x = '', y = 'Genetic variance') +
  facet_wrap( ~ paste(n0, low.var, sep = ', '), ncol = 4) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.background = element_rect(colour = 'black'),
        legend.position = 'none')

fixation.pos.plot = ext.fix %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = p.pos,
      group = interaction(n.pop0, low.var, alpha, extinct),
      linetype = extinct,
      colour = alpha
    ),
    size = 1.5
  ) +
  geom_ribbon(
    aes(
      ymin = p.pos - 2 * sqrt(v.pos / n),
      ymax = p.pos + 2 * sqrt(v.pos / n),
      group = interaction(n.pop0, low.var, alpha, extinct),
      fill = alpha
    ),
    alpha = 0.2
  ) +
  ylim(c(0, 0.5)) +
  scale_colour_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap( ~ paste(n0, low.var, sep = ', '), ncol = 4) +
  guides(linetype = guide_legend('', nrow = 1),
         colour = guide_legend('', nrow = 1),
         fill = guide_legend('', nrow = 1)) +
  labs(x = 'Generation', y = 'Loci at positive fixation') +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = 'none',
        strip.background = element_rect(colour = 'black'),
        strip.text = element_text(size = 12))

fixation.neg.plot = ext.fix %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = p.neg,
      group = interaction(n.pop0, low.var, alpha, extinct),
      linetype = extinct,
      colour = alpha
    ),
    size = 1.5
  ) +
  geom_ribbon(
    aes(
      ymin = p.neg - 2 * sqrt(v.neg / n),
      ymax = p.neg + 2 * sqrt(v.neg / n),
      group = interaction(n.pop0, low.var, alpha, extinct),
      fill = alpha
    ),
    alpha = 0.2
  ) +
  ylim(c(0, 0.5)) +
  scale_colour_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap( ~ paste(n0, low.var, sep = ', '), ncol = 4) +
  guides(linetype = guide_legend('', nrow = 1),
         colour = guide_legend('', nrow = 1),
         fill = guide_legend('', nrow = 1)) +
  labs(x = 'Generation', y = 'Loci at negative fixation') +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = 'none',
        strip.background = element_rect(colour = 'black'),
        strip.text = element_text(size = 12))


varn.legend = get_legend(
  variance.plot + 
    guides(fill = guide_legend(''), colour = guide_legend(''), linetype = guide_legend('')) +
    theme(legend.position = 'bottom',
          legend.box.margin = margin(1, 0, 0, 0))
)

save_plot(plot_grid(variance.plot, fixation.pos.plot, fixation.neg.plot, 
                    varn.legend, rel_heights = c(1, 1, 1, 0.1), ncol = 1), 
          #varn.legend,
          #rel_heights = c(1, 0.1),
          filename = 'simulations/analysis_results/figure_drafts/draft_figs/fig_4.pdf',
          base_width = 8, base_height = 6)

# Summary stats of above:

# Fold-change in fixation of both alleles
ext.fix %>%
  filter(gen %in% 15) %>%
  mutate(ndd = ifelse(alpha %in% 'Density dependent', 'ndd', 'di')) %>%
  select(-c(v.pos, v.neg, n, alpha)) %>%
  pivot_wider(names_from = ndd, values_from = c(p.pos, p.neg)) %>%
  mutate(fold.pos = p.pos_ndd / p.pos_di,
         fold.neg = p.neg_ndd / p.neg_di)

# Change in number of loci changed
ext.fix %>%
  filter(gen %in% 15) %>%
  mutate(ndd = ifelse(alpha %in% 'Density dependent', 'ndd', 'di')) %>%
  select(-c(v.pos, v.neg, n, alpha)) %>%
  pivot_wider(names_from = ndd, values_from = c(p.pos, p.neg)) %>%
  mutate(fold.pos = 25 * (p.pos_ndd - p.pos_di),
         fold.neg = 25 * (p.neg_ndd - p.neg_di))


# Below: old and obsolete plots (code won't work)

# Fixation of only extant populations

ext.fix %>%
  filter(!extinct) %>%
  gather(key = fixtype, value = p, p.fix.pos, p.fix.neg) %>%
  mutate(fixtype = factor(fixtype, labels = c('Negative allele', 'Positive allele'))) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = p,
      group = interaction(n.pop0, low.var, alpha, fixtype),
      linetype = paste(n0, low.var, sep = ', '),
      colour = alpha
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
        strip.text = element_text(size = 12)) #+
  ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_5_extant.pdf',
         height = 5, width = 8)

all.v %>%
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

