# State variables, conditioned on survival/generation of extinction
# (Supplemental figures - "fishscale plots")
# SN - 15 Nov 2021
# (re-run 2 May 2022 lol)

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Clear namespace
rm(list = ls())

# Load in data
all.data = read.csv('simulations/outputs/final_results/alldata_combined.csv')

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
  rename(g = gbar, z = zbar, w = wbar) %>%
  summarise(
    # Population size
    nbar = mean(n),
    nvar = var(n),
    # Genotype
    gbar = mean(g),
    gvar = var(g),
    # Phenotype
    zbar = mean(z),
    zvar = var(z),
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
  labs(x = "Generation", y = "Mean population size") +
  facet_wrap(dd ~ paste(n0, lv, sep = ', '), ncol = 4) +
  theme(
    legend.position = 'bottom',
    panel.grid.major = element_line(colour = 'gray88'),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour = 'black')
    )

ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_supp_ext_gen_n.png',
       width = 8, height = 5)

### Genotype (g)

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
  labs(x = "Generation", y = "Mean population genotype") +
  facet_wrap(dd ~ paste(n0, lv, sep = ', '), ncol = 4) +
  theme(
    legend.position = 'bottom',
    panel.grid.major = element_line(colour = 'gray88'),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour = 'black')
  )

ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_supp_ext_gen_g.png',
       width = 8, height = 5)

### Phenotype (z)

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
      y = zbar,
      group = interaction(ext.gen, extinct),
      colour = extinct
    )
  ) +
  scale_fill_manual(values = c('black', 'red'), '') +
  scale_color_manual(values = c('black', 'red'), '') +
  labs(x = "Generation", y = "Mean population phenotype") +
  facet_wrap(dd ~ paste(n0, lv, sep = ', '), ncol = 4) +
  theme(
    legend.position = 'bottom',
    panel.grid.major = element_line(colour = 'gray88'),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour = 'black')
  )

ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_supp_ext_gen_z.png',
       width = 8, height = 5)

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
  labs(x = "Generation", y = "Mean intrinsic fitness") +
  facet_wrap(dd ~ paste(n0, lv, sep = ', '), ncol = 4) +
  theme(
    legend.position = 'bottom',
    panel.grid.major = element_line(colour = 'gray88'),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour = 'black')
  )

ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_supp_ext_gen_w.png',
       width = 8, height = 5)

### Additive genetic variance

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
  labs(x = "Generation", y = "Mean additive genetic variance") +
  facet_wrap(dd ~ paste(n0, lv, sep = ', '), ncol = 4) +
  theme(
    legend.position = 'bottom',
    panel.grid.major = element_line(colour = 'gray88'),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour = 'black')
  )

ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_supp_ext_gen_v.png',
       width = 8, height = 5)

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
  labs(x = "Generation", y = "Proportion of loci at fixation (positive allele)") +
  facet_wrap(dd ~ paste(n0, lv, sep = ', '), ncol = 4) +
  theme(
    legend.position = 'bottom',
    panel.grid.major = element_line(colour = 'gray88'),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour = 'black')
  )

ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_supp_ext_gen_fixp.png',
       width = 8, height = 5)

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
  labs(x = "Generation", y = "Proportion of loci at fixation (negative allele)") +
  facet_wrap(dd ~ paste(n0, lv, sep = ', '), ncol = 4) +
  theme(
    legend.position = 'bottom',
    panel.grid.major = element_line(colour = 'gray88'),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour = 'black')
  )

ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_supp_ext_gen_fixn.png',
       width = 8, height = 5)
