# Figure 4, showing mean additive genetic variance and fixation rates for 15
# generation simulations (conditioned on extinction/survival to end of sim)
# NOTE: this script also contains code for Fig SI1, showing rates of fixation

# Clear namespace
rm(list = ls())

# Load packages
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

# Load data for intrinsic fitness (W)

### Read in and aggregate data

# All data
all.data = read.csv('simulations/outputs/alldata_combined.csv') %>%
  # Re-start time at generation 0 (instead of generation 1)
  mutate(gen = gen - 1)

# Conditioning genetic variance (v) by extinct/extant
ext.v = all.data %>%
  group_by(n.pop0, low.var, alpha, gen, extinct) %>%
  summarise(vbar = mean(v),
            vvar = var(v),
            n = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("High diversity", "Low diversity")))

# Fixation probabilities by extinction
ext.fix = all.data %>%
  group_by(n.pop0, low.var, alpha, gen, extinct) %>%
  summarise(p.pos = mean(25*p.fix.pos),
            v.pos = var(25*p.fix.pos),
            p.neg = mean(25*p.fix.neg),
            v.neg = var(25*p.fix.neg),
            n = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("High diversity", "Low diversity")))

### Combined figure 4 where panel (a) is genetic variation, (b) and (c) are fixation rates

# Part (a)
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
    linewidth = 1.25
  ) +
  geom_ribbon(
    aes(
      ymin = vbar - 2 * sqrt(vvar / n),
      ymax = vbar + 2 * sqrt(vvar / n),
      group = interaction(n.pop0, low.var, alpha, extinct),
      fill = alpha
    ),
    alpha = 0.2
  ) +
  scale_y_continuous(limits = c(0, 0.5), breaks = (0:5)/10) +
  scale_colour_manual(values = c('black', 'purple'), '') +
  scale_fill_manual(values = c('black', 'purple'), '') +
  scale_linetype_manual(values = c(2, 1), '') +
  labs(x = 'Generation', 
       y = expression(paste('Genetic variance, ', group(langle, (sigma[a]^2)[t], rangle)))) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) +
  facet_wrap( ~ paste(n0, low.var, sep = ', '), ncol = 4) +
  theme(panel.background = element_rect(fill = 'white'),
        # panel.grid = element_line('gray88'),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(colour = 'black'),
        legend.position = 'top')

variance.plot

ggsave('analysis_results/figures/fig4_varloss.png',
       width = 8, height = 3)

### Fixation plots

# Panel (a) - positive fixation
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
  scale_y_continuous(limits = c(0, 13), breaks = (0:4)*3) +
  scale_colour_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_linetype_manual(values = c(1, 2)) +
  facet_wrap( ~ paste(n0, low.var, sep = ', '), ncol = 4) +
  guides(linetype = guide_legend('', nrow = 1),
         colour = guide_legend('', nrow = 1),
         fill = guide_legend('', nrow = 1)) +
  labs(x = '', y = 'Loci at positive fixation') +
  guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.grid.minor = element_blank(),
        legend.position = 'none')
# strip.background = element_rect(colour = 'black'),
# strip.text = element_text(size = 12))

# Part (b) - negative fixation
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
  scale_y_continuous(limits = c(0, 13), breaks = (0:4)*3) +
  scale_colour_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  scale_linetype_manual(values = c(1, 2)) +
  facet_wrap( ~ paste(n0, low.var, sep = ', '), ncol = 4) +
  guides(linetype = guide_legend('', nrow = 1),
         colour = guide_legend('', nrow = 1),
         fill = guide_legend('', nrow = 1)) +
  labs(x = 'Generation', y = 'Loci at negative fixation') +
  guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        strip.text = element_blank(),
        strip.background = element_blank())

# Create legend
varn.legend = get_legend(
  variance.plot + 
    guides(fill = guide_legend(''), colour = guide_legend(''),
           linetype = guide_legend(' ', override.aes = list(linewidth = 0.5), values = c(1, 3))) +
    theme(legend.position = 'bottom',
          legend.box.margin = margin(1, 0, 0, 0))
)

# Fig SI1: fixations of either allele, broken out by extinct/surviving plot

save_plot(plot_grid(fixation.pos.plot, fixation.neg.plot, 
                    labels = c('(A)', '(B)'),
                    label_x = c(0.0125, 0.0125),
                    label_y = c(0.99, 1.02),
                    label_size = 12,
                    varn.legend, rel_heights = c(1, 1, 0.1), ncol = 1),
          filename = 'analysis_results/figures/fig_i1_fixations.png',
          base_width = 8, base_height = 4)

# For combined plot with all three in the same:
# (old version from early version of MS)
save_plot(plot_grid(variance.plot, fixation.pos.plot, fixation.neg.plot, 
                    labels = c('(A)', '(B)', '(C)'),
                    label_x = c(-0.015, -0.015, -0.015),
                    label_y = c(1, 1.16, 1.16),
                    label_size = 12,
                    varn.legend, rel_heights = c(1.2, 1, 1, 0.1), ncol = 1),
          filename = 'analysis_results/figures/fig_4_orig_var_fix.png',
          base_width = 8, base_height = 6)
