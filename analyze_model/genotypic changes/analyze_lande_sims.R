#' ---
#' title: "Comparing Lande with simulations"
#' author: Scott
#' output: github_document
#' ---

#' I have a simulation that should, with a large number of loci, approach the
#' dynamics of phenotypic change described in Lande, 1976. Here, I am testing
#' simulation runs against Lande's expression for several values of $m$ (i.e.,
#' number of loci) to see how the model compares with predictions.
#' 
#' I ran simulation trials at only the very first step in the simulation. This
#' precludes any extinction issues which may bias estimates. I ran with m = 1,
#' 9, 25, 49, and 81 loci and with theta (the initial degree of maladaptation)
#' equal to 2; genotypes should approach normality when the square root of m is
#' much larger than theta.


library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

all.landes = read.csv('analyze_model/genotypic changes/first_step_lande.csv')

lande.summ = all.landes %>%
  group_by(trial, n_loci) %>%
  summarise(d_1 = d_t[1], d_2 = d_t[2],
            lande_d = lande_d_t[1],
            land2_d = land2_d_t[1],
            n_2 = n_t[2],
            v_1 = v_t[1],
            vz1 = vz_t[1])

head(lande.summ)

lande.summ %>%
  filter(n_loci %in% 1) %>%
  ggplot(aes(x = lande_d, y = d_2)) +
  geom_point(aes(colour = n_2)) +
  stat_smooth(method = 'lm', colour = 'purple') +
  geom_segment(aes(x = min(lande_d), xend = max(lande_d),
                   y = min(lande_d), yend = max(lande_d)),
               colour = 'blue') +
  scale_colour_gradient(low = 'gray', high = 'forestgreen') +
  theme(panel.background = element_blank())

# For Brett: all panels, Lande with additive variance
lande.summ %>%
  ggplot(aes(x = lande_d, y = d_2)) +
  geom_point(aes(colour = n_2)) +
  stat_smooth(method = 'lm', colour = 'purple') +
  geom_segment(aes(x = min(lande_d), xend = max(lande_d),
                   y = min(lande_d), yend = max(lande_d)),
               colour = 'blue') +
  scale_colour_gradient(low = 'gray', high = 'forestgreen') +
  theme(panel.background = element_blank()) +
  facet_wrap(~ n_loci)

# Lande expression with phenotypic variance (vzt)
lande.summ %>%
  filter(n_loci %in% 1) %>%
  ggplot(aes(x = land2_d, y = d_2)) +
  geom_point(aes(colour = n_2)) +
  stat_smooth(method = 'lm', colour = 'purple') +
  geom_segment(aes(x = min(lande_d), xend = max(lande_d),
                   y = min(lande_d), yend = max(lande_d)),
               colour = 'blue') +
  scale_colour_gradient(low = 'gray', high = 'forestgreen') +
  theme(panel.background = element_blank())

# For Brett:
lande.summ %>%
  ggplot(aes(x = land2_d, y = d_2)) +
  geom_point(aes(colour = n_2)) +
  stat_smooth(method = 'lm', colour = 'purple') +
  geom_segment(aes(x = min(lande_d), xend = max(lande_d),
                   y = min(lande_d), yend = max(lande_d)),
               colour = 'blue') +
  scale_colour_gradient(low = 'gray', high = 'forestgreen') +
  theme(panel.background = element_blank()) +
  facet_wrap(~ n_loci)

lande.summ %>%
  filter(n_loci > 1) %>%
  ggplot(aes(x = lande_d, y = d_2)) +
  geom_point(aes(colour = v_1)) +
  geom_segment(aes(x = min(lande_d), xend = max(lande_d),
                   y = min(lande_d), yend = max(lande_d)),
               colour = 'blue') +
  stat_smooth(method = 'lm', colour = 'purple') +
  facet_wrap(~ factor(n_loci))

lande.summ %>%
  filter(n_loci > 1) %>%
  ggplot(aes(x = land2_d, y = d_2)) +
  geom_point(aes(colour = v_1)) +
  geom_segment(aes(x = min(lande_d), xend = max(lande_d),
                   y = min(lande_d), yend = max(lande_d)),
               colour = 'blue') +
  stat_smooth(method = 'lm', colour = 'purple') +
  facet_wrap(~ factor(n_loci))

lande.summ %>%
  ggplot(aes(x = v_1, y = vz1)) +
  geom_point(aes(fill = factor(n_loci)),
             shape = 21) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 5, 'BuGn'))
# So... you get much more phenotypic variance than additive genetic variance?

lande.summ %>%
  ggplot(aes(x = lande_d, y = land2_d)) +
  geom_point(aes(colour = sqrt(n_loci)),
             alpha = 0.5) +
  geom_segment(aes(x = 1.2, xend = 2.2,
                   y = 1.2, yend = 2.2),
               colour = 'red')

lande.summ %>%
  ggplot(aes(x = d_2 - lande_d)) +
  geom_histogram(binwidth = 0.05) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 200),
               colour = 'red') +
  facet_wrap(~ n_loci)

lande.summ %>%
  ggplot(aes(x = d_2 - land2_d)) +
  geom_histogram(binwidth = 0.05) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 200),
               colour = 'red') +
  facet_wrap(~ n_loci)

# Okay, so consistently biased high, doesn't seem much better for larger m.

lande.summ %>%
  group_by(n_loci) %>%
  summarise(diff_bar = mean(d_2 - lande_d),
            diff_se = var(d_2 - lande_d) / sqrt(2000))

# No evidence of pattern wth $m$. But yes clearly biased high,
# i.e., Lande underpredicts (why?)

lm(d_2 ~ lande_d * factor(n_loci), lande.summ) %>% summary()
lm(d_2 ~ lande_d * sqrt(n_loci), lande.summ) %>% summary()
lm(d_2 ~ lande_d + sqrt(n_loci), lande.summ) %>% summary()
lm(d_2 ~ lande_d, lande.summ) %>% summary()

lm(d_2 ~ lande_d - 1, lande.summ) %>% summary()
lm(d_2 ~ lande_d - 1, lande.summ) %>% plot()

