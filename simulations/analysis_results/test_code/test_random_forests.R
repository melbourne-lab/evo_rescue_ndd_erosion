# Crude testing of random forests, looking at feature importance
# SN - Sept 7 2020

library(randomForest)

### Running with 4000 trials with all data stored

small.lovar.did = read.csv('simulations/outputs/test_alldata_n20_a000_lowvar.csv')

head(small.lovar.did)

sli.g5 = small.lovar.did %>% 
  filter(ext.gen > 5 & gen < 6) %>%
  select(-c(ext.gen, n.pop0)) %>%
  pivot_wider(
    names_from = gen,
    values_from = c(n, gbar, zbar, wbar, v, p.fix.pos, p.fix.neg),
    values_fn = list(v = mean, p.fix.pos = mean, p.fix.neg = mean,
                     gbar = mean, zbar = mean, wbar = mean, n = mean)
  ) %>%
  mutate(extinct = factor(extinct)) %>%
  select(-c(p.fix.pos_1, p.fix.neg_1, n_1))

set.seed(24924)
test.trials = sample(sli.g5$trial)

sli5.train = sli.g5 %>% filter(trial %in% test.trials[1:(.75*length(test.trials))]) %>% select(-trial)
sli5.test = sli.g5 %>% filter(!trial %in% test.trials[1:(.75*length(test.trials))]) %>% select(-trial)

sli5.1 = randomForest(extinct ~ ., data = sli5.train, importance = TRUE)

sli5.1

sli5.1$importance
plot(sli5.1$importance)
varImpPlot(sli5.1)

sli5.1.valid = predict(sli5.1, sli5.test)
table(sli5.test$extinct, sli5.1.valid) %>% (function(x) (x[1] + x[4]) / sum(x))

### Try with three generations

sli.g3 = small.lovar.did %>% 
  filter(ext.gen > 3 & gen < 4) %>%
  select(-c(ext.gen, n.pop0)) %>%
  pivot_wider(
    names_from = gen,
    values_from = c(n, gbar, zbar, wbar, v, p.fix.pos, p.fix.neg),
    values_fn = list(v = mean, p.fix.pos = mean, p.fix.neg = mean,
                     gbar = mean, zbar = mean, wbar = mean, n = mean)
  ) %>%
  mutate(extinct = factor(extinct)) %>%
  select(-c(p.fix.pos_1, p.fix.neg_1, n_1))

set.seed(133094)
test.trials3 = sample(sli.g3$trial)

sli3.train = sli.g3 %>% filter(trial %in% test.trials[1:(.75*length(test.trials3))]) %>% select(-trial)
sli3.test = sli.g3 %>% filter(!trial %in% test.trials[1:(.75*length(test.trials3))]) %>% select(-trial)

sli3.1 = randomForest(extinct ~ ., data = sli3.train, importance = TRUE)

sli3.1

sli3.1$importance
plot(sli3.1$importance)
varImpPlot(sli3.1)

sli3.1.valid = predict(sli3.1, sli3.test)
table(sli3.test$extinct, sli3.1.valid) %>% (function(x) (x[1] + x[4]) / sum(x))

### Try first three generations for populations ext.gen > 5

sli.g3.5 = small.lovar.did %>% 
  filter(ext.gen > 5 & gen < 4) %>%
  select(-c(ext.gen, n.pop0)) %>%
  pivot_wider(
    names_from = gen,
    values_from = c(n, gbar, zbar, wbar, v, p.fix.pos, p.fix.neg),
    values_fn = list(v = mean, p.fix.pos = mean, p.fix.neg = mean,
                     gbar = mean, zbar = mean, wbar = mean, n = mean)
  ) %>%
  mutate(extinct = factor(extinct)) %>%
  select(-c(p.fix.pos_1, p.fix.neg_1, n_1))

set.seed(2323)
test.trials3.5 = sample(sli.g3.5$trial)

sli3.5.train = sli.g3.5 %>% filter(trial %in% test.trials3.5[1:(.75*length(test.trials3.5))]) %>% select(-trial)
sli3.5.test = sli.g3.5 %>% filter(!trial %in% test.trials3.5[1:(.75*length(test.trials3.5))]) %>% select(-trial)

sli3.5.1 = randomForest(extinct ~ ., data = sli3.5.train, importance = TRUE)

sli3.5.1

sli3.5.1$importance
plot(sli3.5.1$importance)
varImpPlot(sli3.5.1)

sli3.5.1.valid = predict(sli3.5.1, sli3.5.test)
table(sli3.5.test$extinct, sli3.5.1.valid) %>% (function(x) (x[1] + x[4]) / sum(x))

### Old, less refined stuff.
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

### Thoughts on random forests

# Not being able to handle missing data is a massive limitation.
# It means that you need to do censoring - extinction time beyond a certain point.
# So you need to do Pr(extinct after certain time) with data from before that time.
# Seems difficult.
# 
# It looks like g_{t-1} has the best predictive ability (highest decrease in Gini).
# More so than population size (n_{t-1})
# Fixation probabilities have close to no predictive power
# Although, during the first several time steps, there are no fixations.
#
# Overall, I'm not optimistic about prospects here.
# This also seems to be true of any linear modeling based approaches -
# because those also can't handle missing data.


