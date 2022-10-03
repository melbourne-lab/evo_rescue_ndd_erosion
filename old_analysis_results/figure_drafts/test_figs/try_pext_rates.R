# Script to try to get a P(extinct) plots against rates

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

rm(list = ls())

all.data = read.csv('simulations/outputs/final_results/alldata_combined.csv')

all.tau = all.data %>%
  group_by(trial, n.pop0, low.var, alpha) %>%
  mutate(loglam = c(NA, diff(log(n))),
         d = abs(2.75 - gbar),
         k = c(NA, exp(diff(log(d)))),
         varred = c(NA, exp(diff(log(v))))) %>%
  mutate(tau = max(gen) - gen + 1 - as.numeric(any(n < 2)),
         ext = max(gen) < 15 | any(n < 2),
         ext5 = tau < 5) %>%
  filter(!is.na(loglam), gen > 1) %>%
  ungroup() %>%
  mutate(n.pop0 = factor(n.pop0),
         alpha  = factor(alpha))

ll.blank = expand.grid(
  loglam = (-7:4) / 10,
  low.var = c(TRUE, FALSE),
  n.pop0 = factor(c(20, 100)),
  alpha = factor(c(0, 0.0035))
)

llmods = all.tau %>%
  filter(gen < 11) %>%
  split(.$gen) %>%
  map(~ glm(ext5 ~ low.var * n.pop0 * alpha * loglam,
            data = .,
            family = 'binomial'))

ll.preds = lapply(llmods, function(x) predict.glm(object = x, newdata = ll.blank, type = 'resp')) %>%
  do.call(what = cbind) %>%
  cbind(ll.blank, .) %>%
  gather(key = gen, value = p.ext, -c(loglam, low.var, n.pop0, alpha)) %>%
  mutate(gen = as.numeric(gen))
  
ggplot(ll.preds) +
  geom_tile(
    aes(
      x = gen,
      y = loglam,
      fill = p.ext
    )
  ) +
  facet_wrap(n.pop0 + low.var ~ alpha, nrow = 4) +
  scale_fill_viridis_c()

ll.r2 = llmods %>%
  lapply(function(x) x %>%
           predict.glm(type = 'resp') %>%
           cbind(x$data, p.ext = .)) %>%
  do.call(what = rbind) %>%
  mutate(pred.ext = p.ext > 0.5) %>%
  group_by(gen, n.pop0, low.var, alpha) %>%
  summarise(prob.correct = mean(pred.ext == ext5))

ll.plot = merge(ll.preds, ll.r2)

ll.plot %>% 
  filter(gen < 15) %>%
  mutate(low.var = factor(low.var, labels = c('high div.', 'low div.')),
         n.pop0 = factor(n.pop0, labels = c('small', 'large')),
         alpha = factor(alpha, labels = c('dens. indep.', 'dens. dep.'))) %>%
  ggplot() +
  geom_tile(
    aes(
      x = gen,
      y = loglam,
      fill = p.ext,
      alpha = prob.correct
    )
  ) +
  facet_wrap(paste(n.pop0, low.var, sep = ', ') ~ alpha, nrow = 4) +
  scale_alpha_continuous(range = c(min(ll.r2$prob.correct), 1)) +
  scale_fill_viridis_b(option = 'A') +
  theme(panel.grid = element_blank(),
        panel.background = element_blank())

### Try a plot that compares NDD and DI in one panel

ll.linpreds = lapply(llmods, function(x) predict.glm(object = x, newdata = ll.blank)) %>%
  do.call(what = cbind) %>%
  cbind(ll.blank, .) %>%
  gather(key = gen, value = p.ext, -c(loglam, low.var, n.pop0, alpha)) %>%
  mutate(gen = as.numeric(gen)) %>%
  spread(key = alpha, value = p.ext) %>%
  mutate(d.alpha = `0.0035` - `0`) %>%
  select(-c(`0`, `0.0035`))

ggplot(ll.linpreds) +
  geom_tile(
    aes(
      x = gen,
      y = loglam,
      fill = d.alpha
    )
  ) +
  facet_wrap(~ n.pop0 + low.var, nrow = 4) +
  scale_fill_viridis_c()

# Scale here is whack.
