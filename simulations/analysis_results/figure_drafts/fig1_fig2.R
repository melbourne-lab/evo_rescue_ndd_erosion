### Population size and dynamics figures.

# Load packages

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

# Load in all simulation data
all.data = read.csv('simulations/outputs/final_results/alldata_combined.csv')

##### Figure 1


### Aggregate mean population size (and variance for standard error)

all.n = all.data %>%
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
  # Get means and variances
  group_by(gen, n.pop0, low.var, alpha) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            n.trials = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Initially small", "Initially large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("Low genetic diversity", "High genetic diversity")))

all.n %>%
  mutate(n0 = factor(n0, levels = levels(n0)[2:1])) %>%
  ggplot(aes(x = gen)) +
  geom_segment(
    aes(
      x = 1, xend = 15,
      y = n.pop0, yend = n.pop0
    ),
    linetype = 3,
    size = 0.5,
    colour = 'gray'
  ) +
  geom_line(
    aes(
      y = nbar,
      group = interaction(alpha, n.pop0, low.var),
      linetype = low.var,
      colour = alpha
    ),
    size = 1.25
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / n.trials),
      ymax = nbar + 2 * sqrt(nvar / n.trials),
      group = interaction(alpha, n.pop0, low.var),
      fill = alpha
    ),
    alpha = 0.2,
    size = 0.125
  ) +
  scale_color_manual(
    values = c('black', 'purple'),
    labels = c("Density\nindependent", "Density\ndependent"),
    name = ""
  ) +
  scale_fill_manual(
    values = c('black', 'purple'),
    labels = c("Density\nindependent", "Density\ndependent"),
    name = ""
  ) +
  scale_linetype(# _manual(
    # values = c(1, 5),
    name = "",
    labels = c("High diversity", "Low diversity")
  ) +
  labs(x = 'Generation', y = 'Mean population size') +
  scale_y_log10() +
  facet_wrap(~ n0) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.direction = 'horizontal',
        legend.position = c(0.2, 0.85),
        legend.text = element_text(size = 12),
        strip.background = element_rect(colour = 'black'),
        strip.text = element_text(size = 12)) #+

  ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_pop_size_large_first.png',
         width = 8, height = 5)

##### Figure S1

### Get mean population size by time to extinction (or extant)

ext.n = all.data %>%
  group_by(gen, n.pop0, low.var, alpha, ext.gen) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            n.trials = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Initially small", "Initially large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("High genetic diversity", "Low genetic diversity")))

### Plot

ext.n %>%
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
    data = . %>% filter(ext.gen < 15),
    aes(
      y = nbar,
      group = interaction(alpha, n.pop0, low.var, ext.gen)
    ),
    colour = 'firebrick3',
    size = 0.5
  ) +
  geom_line(
    data = . %>% filter(ext.gen == 15),
    aes(
      y = nbar,
      group = interaction(alpha, n.pop0, low.var)
    ),
    colour = 'black',
    size = 1
  ) +
  geom_ribbon(
    data = . %>% filter(ext.gen < 15),
    aes(
      ymin = nbar - sqrt(nvar / n.trials),
      ymax = nbar + sqrt(nvar / n.trials),
      group = interaction(alpha, n.pop0, low.var, ext.gen)
    ),
    fill = 'firebrick3',
    alpha = 0.2,
    size = 0.125
  ) +
  geom_ribbon(
    data = . %>% filter(ext.gen == 15),
    aes(
      ymin = nbar - sqrt(nvar / n.trials),
      ymax = nbar + sqrt(nvar / n.trials),
      group = interaction(alpha, n.pop0, low.var)
    ),
    fill = 'black',
    alpha = 0.2,
    size = 0.125
  ) +  
  labs(x = 'Generation', y = 'Mean population size') +
  scale_y_log10() +
  facet_grid(cols = vars(reorder(n0, desc(n0)), low.var), rows = vars(alpha)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
#        legend.direction = 'horizontal',
#        legend.position = c(0.2, 0.8),
        legend.text = element_text(size = 12),
        strip.background = element_rect(colour = 'black'),
        strip.text = element_text(size = 12)) #+
  ggsave('simulations/analysis_results/figure_drafts/draft_figs/fig_s1.pdf',
         width = 8, height = 5)

##### Figure 2

all.ext = all.data %>%
  # Get number of extant trials in each generation for each treatment
  group_by(n.pop0, low.var, alpha, gen) %>%
  summarise(n.trials = n()) %>%
  # Summarise extinctions
  mutate(# Number of trials going extinct
         d.trials = -c(0, diff(n.trials)),
         # Proportion of extant trials going extinct
         prop.extant = c(d.trials[-1], 0) / n.trials,
         # Proportion of all trials going extinct
         prop.totaln = c(d.trials[-1], 0) / max(n.trials),
         # Proportion of extinct trials extinct in this generation
         prop.condit = c(d.trials[-1], 0) / sum(d.trials),
         # Proportion of trials surviving
         prop.surviv = n.trials / max(n.trials)) %>%
  ungroup() %>%
  mutate(n.pop0 = factor(n.pop0, labels = c("Small", "Large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("Low diversity", "High diversity")))

### Get genotype vs. extinction risk

ext.probs = all.data %>% distinct(n.pop0, alpha, low.var, trial, .keep_all = TRUE)

# fit model
all.ext.glm.3 = glm(extinct ~ factor(n.pop0) * low.var * factor(alpha) +
                              factor(n.pop0) * factor(alpha) * gbar, 
                    data = ext.probs,
                    family = 'binomial')

predicted.surv = expand.grid(low.var = c(TRUE, FALSE),
                             n.pop0 = factor(c(20, 100)),
                             alpha = factor(c(0, 0.0035)),
                             gbar = (-5:5)/10) %>%
  cbind(logit.extinct = predict(all.ext.glm.3, newdata = .)) %>%
  mutate(p.extinct = logit.extinct %>% (function(x) exp(x) / (1 + exp(x)))) %>%
  mutate(low.var = factor(low.var, labels = c("High diversity", "Low diversity")),
         n.pop0 = factor(n.pop0, labels = c("Small", "Large")))

# Instant probability of extinction
inst.probs = ggplot(all.ext %>% filter(gen < 15), aes(x = gen)) +
  geom_line(
    aes( 
      y = prop.extant, 
      group = interaction(n.pop0, low.var, alpha),
      colour = alpha
    )
  ) +
  geom_ribbon(
    aes(
      ymin = prop.extant - 2 * sqrt(prop.extant * (1-prop.extant) / n.trials),
      ymax = prop.extant + 2 * sqrt(prop.extant * (1-prop.extant) / n.trials),
      group = interaction(n.pop0, low.var, alpha),
      fill = alpha
    ),
    alpha = 0.2
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1) +
  labs(x = 'Generation', y = '') +
  theme(legend.position = 'none',
        panel.background = element_blank(),
        strip.background = element_rect(colour = 'black'))

# Cumulative probability of extinction
cuml.probs = ggplot(all.ext %>% filter(gen < 15), aes(x = gen)) +
  geom_ribbon(
    aes(
      xmin = 0, xmax = 14,
      ymin = 0, ymax = 1 - prop.surviv,
      group = interaction(n.pop0, low.var, alpha),
      fill = factor(alpha)
    ),
    alpha = 0.25
  ) +
  scale_fill_manual(values = c('black', 'purple')) +
  labs(x = 'Generation', y = 'Probability of extinction') +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.background = element_blank(),
        strip.background = element_rect(colour = 'black'))

# Genotype vs. probability of extinction
geno.probs = ggplot(predicted.surv) +
  geom_line(
    aes(
      x = gbar, 
      y = p.extinct, 
      group = interaction(n.pop0, low.var, alpha), 
      colour = alpha
    )
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  labs(x = 'Initial genotype', y = '') +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 1) +
  theme(legend.position = 'none',
        panel.background = element_blank(),
        strip.background = element_rect(colour = 'black'))

extinct.legend = get_legend(
  inst.probs + 
    guides(fill = guide_legend(''), colour = guide_legend('')) +
    theme(legend.position = 'bottom',
          legend.box.margin = margin(1, 0, 0, 0))
)

data.plots = plot_grid(cuml.probs, inst.probs, geno.probs, 
                       labels = c('(A)', '(B)', '(C)'),
                       nrow = 1)
plot_grid(data.plots, extinct.legend, ncol = 1, rel_heights = c(1, .1))

plot_grid(data.plots, extinct.legend, ncol = 1, rel_heights = c(1, .1)) %>%
  save_plot(filename = 'simulations/analysis_results/figure_drafts/draft_figs/fig_2.pdf',
            base_width = 8, base_height = 8)

##### Figure S2

### Get population growth rates of _extant/recently extinct_ populations

all.lambda = all.data %>%
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
  # group_by(trial, n.pop0, low.var, alpha) %>%
  # filter(any(n>0)) %>%
  # Get rid of unnecessary zeros (more than one gen post-extincton)
  # start by adding difference column
  group_by(trial, n.pop0, low.var, alpha) %>%
  mutate(ndiff = c(0, diff(n))) %>%
  # get rid of rows where population is extinct and size did not change
  filter(! (!n & !ndiff)) %>%
  # Get growth rate
  mutate(lambda = c(exp(diff(log(n))), NA)) %>%
  # Get rid of last time step (set to lambda = NA above)
  filter(!is.na(lambda)) %>%
  # means and variances of lambda
  group_by(gen, n.pop0, low.var, alpha) %>%
  summarise(lambda.bar = mean(lambda),
            lambda.var = var(lambda),
            n.trials = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Initially small", "Initially large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("Low genetic variance", "High genetic variance")))

### Extinction plot

all.lambda %>%
  ggplot(aes(x = gen)) +
  geom_segment(
    aes(
      x = 1, xend = 14,
      y = 1, yend = 1
    ),
    linetype = 3,
    colour = 'lightgray'
  ) +
  geom_line(
    aes(
      y = lambda.bar,
      group = interaction(alpha, n.pop0, low.var),
      linetype = low.var,
      colour = alpha
    ),
    size = 1.25
  ) +
  geom_ribbon(
    aes(
      ymin = lambda.bar - 2 * sqrt(lambda.var / n.trials),
      ymax = lambda.bar + 2 * sqrt(lambda.var / n.trials),
      group = interaction(alpha, n.pop0, low.var),
      fill = alpha
    ),
    alpha = 0.2,
    size = 0.125
  ) +
  scale_color_manual(
    values = c('black', 'purple'),
    labels = c("Density\nindependent", "Density\ndependent"),
    name = ""
  ) +
  scale_fill_manual(
    values = c('black', 'purple'),
    labels = c("Density\nindependent", "Density\ndependent"),
    name = ""
  ) +
  scale_linetype(
    name = "",
    labels = c("High variation", "Low variation")
  ) +
  labs(x = 'Generation', y = 'Mean population size') +
  facet_wrap(~ n0) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        # legend.direction = 'horizontal',
        # legend.position = c(0.2, 0.8),
        legend.position = 'bottom',
        legend.text = element_text(size = 12),
        strip.background = element_rect(colour = 'black'),
        strip.text = element_text(size = 12))

# Hm - is this worth including?

### Old code

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

### A plot of lambda (population growth rate)

all.lambda = read.csv('simulations/outputs/final_results/alldata_combined.csv') %>%
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
  # Get rid of unnecessary zeros (more than one gen post-extincton)
  # start by adding difference column
  mutate(ndiff = c(0, diff(n))) %>%
  # get rid of rows where population is extinct and size did not change
  filter(! (!n & !ndiff)) %>%
  # Get growth rate
  mutate(lambda = c(exp(diff(log(n))), NA)) %>%
  # Get rid of last time step (set to lambda = NA above)
  filter(!is.na(lambda)) %>%
  # means and variances of lambda
  group_by(gen, n.pop0, low.var, alpha) %>%
  summarise(lambda.bar = mean(lambda),
            lambda.var = var(lambda),
            n.trials = n()) %>%
  ungroup() %>%
  mutate(n0 = factor(n.pop0, labels = c("Initially small", "Initially large")),
         alpha = factor(alpha, labels = c("Density independent", "Density dependent")),
         low.var = factor(low.var, labels = c("Low genetic variance", "High genetic variance")))

### Lambda plot

all.lambda %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = lambda.bar,
      group = interaction(alpha, n.pop0, low.var),
      linetype = low.var,
      colour = alpha
    ),
    size = 1.25
  ) +
  geom_ribbon(
    aes(
      ymin = lambda.bar - 2 * sqrt(lambda.var / n.trials),
      ymax = lambda.bar + 2 * sqrt(lambda.var / n.trials),
      group = interaction(alpha, n.pop0, low.var),
      fill = alpha
    ),
    alpha = 0.2,
    size = 0.125
  ) +
  scale_color_manual(
    values = c('black', 'purple'),
    labels = c("Density\nindependent", "Density\ndependent"),
    name = ""
  ) +
  scale_fill_manual(
    values = c('black', 'purple'),
    labels = c("Density\nindependent", "Density\ndependent"),
    name = ""
  ) +
  scale_linetype(
    name = "",
    labels = c("High variation", "Low variation")
  ) +
  labs(x = 'Generation', y = 'Mean population size') +
  facet_wrap(~ n0) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        # legend.direction = 'horizontal',
        # legend.position = c(0.2, 0.8),
        legend.position = 'bottom',
        legend.text = element_text(size = 12),
        strip.background = element_rect(colour = 'black'),
        strip.text = element_text(size = 12))

### Some summary stats

# Size at end of simulation
all.n %>%
  filter(gen %in% 15) %>%
  mutate(ndd = ifelse(alpha %in% 'Density dependent', 'ndd', 'di')) %>%
  select(-c(alpha, nvar, n.trials, n0)) %>%
  spread(key = ndd, value = nbar) %>%
  mutate(ratio.n = 1 - (ndd / di),
         ldiff.n = log(di) - log(ndd))
  
all.n.nosum = read.csv('simulations/outputs/final_results/alldata_combined.csv') %>%
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
  filter(any(n>0)) 

all.n.nosum %>%
  filter(gen %in% 15) %>%
  glm(formula = n ~ factor(n.pop0) * factor(alpha) * low.var,
      family = 'poisson') %>%
  summary()

# At minimum size

all.n %>%
  group_by(n.pop0, alpha, low.var) %>%
  filter(nbar %in% min(nbar)) %>%
  ungroup() %>%
  mutate(ndd = ifelse(alpha %in% 'Density dependent', 'ndd', 'di')) %>%
  select(-c(alpha, nvar, n.trials, n0, gen)) %>%
  spread(key = ndd, value = nbar) %>%
  mutate(ratio.n = 1 - (ndd / di),
         ldiff.n = log(di) - log(ndd))
