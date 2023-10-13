# Script to generate robustness plots for environmental contributions to
# phenotypic variance (sigma^2_e)

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

rm(list = ls())

# Read in file
sim.out = read.csv('simulations/outputs/robust_checks_out/sig2e_robust_check.csv')

# Check to see it looks correct
head(sim.out)

# Adding rows which treat extinctions as zeros
sim.out = sim.out %>%
  mutate(gen = paste0('gen', gen)) %>%
  pivot_wider(names_from = gen, values_from = n, values_fill = 0) %>%
  pivot_longer(cols = starts_with('gen'), names_to = 'gen', values_to = 'n') %>%
  mutate(gen = gsub('gen', '', gen) %>% as.numeric()) %>%
  arrange(trial, gen)

# Checking to see if any zeros were added inappropriately (size-censoring)
sim.out %>%
  group_by(trial) %>%
  filter(any(!n)) %>%
  mutate(nd = c(NA, diff(n))) %>%
  filter(nd < -1000)
# great - nothing cut out prematurely!

# How many rows (should be a neat even number)
nrow(sim.out)

# Get a summary data frame
# outputs: mean pop size, variance in pop size, survival/extinciton rates, sample size
sim.sum = sim.out %>%
  group_by(gen, n.pop0, lowvar, alpha, sig.e) %>%
  summarise(
    nbar = mean(n), nvar = var(n), 
    pext = mean(!n), nsurv = sum(n > 0), 
    nn = n()
  ) %>%
  ungroup() %>%
  mutate(
    n0 = factor(n.pop0, labels = c("Small", "Large")),
    alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
    lowvar = factor(lowvar, labels = c("High diversity", "Low diversity")),
    sig2e = round(sig.e^2, 2)
  )

head(sim.sum)

# A basic plot of population size
plot.top = sim.sum %>%
  mutate(alpha = factor(alpha)) %>%
  filter(sig2e %in% 0.25) %>%
  ggplot(aes(x = gen)) +
  geom_segment(
    aes(
      x = 0, xend = 15,
      y = n.pop0, yend = n.pop0
    ),
    linetype = 3,
    linewidth = 0.5,
    colour = 'gray'
  ) +
  geom_line(
    aes(
      y = nbar,
      group = alpha,
      colour = alpha
    )
  ) +
  labs(x = '', y = '') +
  facet_grid( ~ paste0(n0, ', ', lowvar)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_log10(
    limits = 10^(.55 + c(0, 1)*2.5),
    sec.axis = sec_axis(~ ., name = expression(sigma[e]^2 == 0.25), labels = NULL, breaks = NULL)
  ) +
  scale_colour_manual(values = c('black', 'purple'), '') +
  scale_fill_manual(values = c('black', 'purple'), '') +
  theme(
    axis.line.y = element_line(),
    axis.line.y.right = element_blank(),
    axis.text.y = element_text(margin = margin(l = 0.09, unit = 'cm')),
    axis.title.y.right = element_text(size = rel(0.8)),
    legend.position = 'none',
    panel.background = element_blank()
  )

plot.mid = sim.sum %>%
  mutate(alpha = factor(alpha)) %>%
  filter(sig2e %in% 1) %>%
  ggplot(aes(x = gen)) +
  geom_segment(
    aes(
      x = 0, xend = 15,
      y = n.pop0, yend = n.pop0
    ),
    linetype = 3,
    linewidth = 0.5,
    colour = 'gray'
  ) +
  geom_line(
    aes(
      y = nbar,
      group = alpha,
      colour = alpha
    )
  ) +
  labs(x = '', y = '') +
  facet_grid(~ paste0(n0, ', ', lowvar)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_log10(
    limits = 10^(0.5 + c(0, 1)*2.5),
    sec.axis = sec_axis(~ ., name = expression(sigma[e]^2 == 1), labels = NULL, breaks = NULL)
  ) +
  scale_colour_manual(values = c('black', 'purple'), '') +
  scale_fill_manual(values = c('black', 'purple'), '') +
  theme(
    axis.line.y = element_line(),
    axis.line.y.right = element_blank(),
    #axis.text.y = element_text(margin = margin(l = 0.18, unit = 'cm')),
    axis.title.y.right = element_text(size = rel(0.8)),
    legend.position = 'none',
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )

# Bottom row
plot.bot = sim.sum %>%
  mutate(alpha = factor(alpha)) %>%
  filter(sig2e %in% 2) %>%
  ggplot(aes(x = gen)) +
  geom_segment(
    aes(
      x = 0, xend = 15,
      y = n.pop0, yend = n.pop0
    ),
    linetype = 3,
    linewidth = 0.5,
    colour = 'gray'
  ) +
  geom_line(
    aes(
      y = nbar,
      group = alpha,
      colour = alpha
    )
  ) +
  labs(x = '', y = '') +
  facet_grid( ~ paste0(n0, ', ', lowvar)) +
  scale_y_log10(
    limits = 10^(.25 + c(0, 1)*2.5),
    breaks = c(10, 100),
    sec.axis = sec_axis(~ ., name = expression(sigma[e]^2 == 2), labels = NULL, breaks = NULL)
  ) +
  scale_colour_manual(values = c('black', 'purple'), '') +
  scale_fill_manual(values = c('black', 'purple'), '') +
  theme(
    axis.line = element_line(),
    axis.line.y.right = element_blank(),
    axis.text.y = element_text(margin = margin(l = 0.225, unit = 'cm')),
    axis.title.y.right = element_text(size = rel(0.8)),
    legend.position = 'none',
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )

plot.leg = get_legend(plot.top + theme(legend.position = 'top'))

n.plot = plot_grid(
  plot.leg, plot.top, plot.mid, plot.bot, 
  nrow = 4, rel_heights = c(.2, 1.15, 1, 1.1)# ,
  # labels = c('', 'a)', 'b)', 'c)'), label_x = -.005
)

gridExtra::grid.arrange(
  gridExtra::arrangeGrob(
    n.plot,
    left = grid::textGrob("Mean population size", rot = 90, vjust = 1.5),
    bottom = grid::textGrob("Generation", vjust = -0.5)
  )
) %>%
  save_plot(
    filename = "analysis_results/robustness_checks/figs/fig_robust_sig2e_popsize.png",
    base_height = 5, base_width = 8
  )

# A plot of extinction rates over time
sim.sum %>%
  # Get proper indexing for extinctions
  filter(gen > 0) %>% mutate(gen = gen - 1) %>%
  # Labels for arranging plot
  mutate(
    alpha = factor(alpha), 
    facet.x = paste0(n0, ', ', lowvar),
    facet.y = paste("sigma[e]^2 ==", sig2e)
  ) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = pext,
      group = alpha,
      colour = alpha
    ),
    linewidth = 0.1
  ) +
  geom_ribbon(
    aes(ymin = 0, ymax = pext, group = alpha, fill = alpha),
    alpha = 0.2
  ) +
  geom_line(
    data = sim.sum %>%
      arrange(desc(gen)) %>% 
      group_by(n.pop0, lowvar, alpha, sig2e) %>% 
      mutate(pins.ext = c(NA, diff(nsurv)) / nsurv) %>% 
      filter(nsurv > 10, !is.na(pins.ext)) %>%
      mutate(
        alpha = factor(alpha), 
        facet.x = paste0(n0, ', ', lowvar),
        facet.y = paste("sigma[e]^2 ==", sig2e)
      ),
    aes(y = pins.ext, group = alpha, colour = alpha)
  ) +
  labs(x = 'Generation', y = 'Probability of extinction') +
  facet_grid(
    facet.y ~ facet.x, 
    labeller = labeller(facet.x = label_value, facet.y = label_parsed)
  ) +
  scale_colour_manual(
    values = c('purple', 'black'), 
    breaks = c("Density dependent", "Density independent"),
    'Instantaneous: '
  ) +
  scale_fill_manual(
    values = c('purple' ,'black'), 
    breaks = c("Density dependent", "Density independent"),
    '  Cumulative: '
  ) +
  guides(
    colour = guide_legend(direction = 'vertical', title.position = 'left'), 
    fill = guide_legend(direction = 'vertical', title.position = 'left')
  ) +
  theme(
    legend.position = 'top',
    legend.title.align = 0.5,
    legend.title = element_text(size = rel(0.8)),
    axis.line = element_line(),
    panel.background = element_blank(),
    strip.text.y = element_text(size = rel(1)),
    strip.background.y = element_blank()
  )

ggsave(
  'analysis_results/robustness_checks/figs/fig_robust_sig2e_extinct.png',
  width = 8, height = 5
)
