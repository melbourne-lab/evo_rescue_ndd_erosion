# Crude testing of random forests, looking at feature importance
# SN - Sept 7 2020

library(randomForest)

popndd = read.csv('simulations/outputs/test_gridded_v.csv') %>% mutate(ndd = TRUE)
popdid = read.csv('simulations/outputs/test_gridded_v2.csv') %>% mutate(ndd = FALSE)

# first try did

popdid.small.g5 = popdid %>% 
  filter(ext.gen > 5 & gen < 6 & n.pop0 %in% 20) %>%
  select(-c(n.pop0, extinct, ndd)) %>%
  pivot_wider(
    names_from = gen,
    values_from = c(v, p.fix.pos, p.fix.neg),
    values_fn = list(v = mean, p.fix.pos = mean, p.fix.neg = mean)
    ) %>%
  mutate(ext.gen = factor(ext.gen)) %>%
  select(-c(p.fix.pos_1, p.fix.neg_1))

set.seed(24924)
test.trials = sample(popdid.small.g5$trial)

popdid.train = popdid.small.g5 %>% filter(trial %in% test.trials[1:3000]) %>% select(-trial)
popdid.test = popdid.small.g5 %>% filter(!trial %in% test.trials[1:3000]) %>% select(-trial)

apply(popdid.train, 2, function(x) sum(is.na(x)))

set.seed(999)
test.small = randomForest(ext.gen ~ ., data = popdid.train)

test.small$importance

# how on earth does v_1 have more importance than p.fig.neg_5??

ggplot(popdid.train) +
  geom_point(aes(x = v_1, y = ext.gen),
             position = position_jitter(height = 0.25))

ggplot(popdid.train) +
  geom_point(aes(x = p.fix.neg_5, y = ext.gen),
             position = position_jitter(width = 0.01, height = 0.25),
             alpha = 0.5)
# huh... interesting! not a very strong effect here

ggplot(popdid.train) +
  geom_point(aes(x = v_5, y = ext.gen),
             position = position_jitter(height = 0.25))
# but also no strong effects here...