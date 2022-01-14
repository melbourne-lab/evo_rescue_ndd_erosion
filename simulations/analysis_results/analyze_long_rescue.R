### Script assessing probability of rescue in long (50-generation) simulations
### SN - 14 Jan 2022

### Setup

# Load packages
library(rstanarm)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

# Clear namespace
rm(list = ls())

# Declare auxiliary function for determining minimum generation at which there
# are three straight generations of some event.

which.diff = function(x, lag) {
  if(length(x) >= lag) { x[which.min(diff(x, lag = (lag - 1)) == (lag-1))] 
  } else { 0 }
}

# Load in simulation data
all.data = 
  rbind(
    read.csv('simulations/outputs/longsims/longsims_n100_a000_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n100_a000_lowvar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a000_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a000_lowvar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n100_a035_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n100_a035_lowvar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a035_hivar.csv'),
    read.csv('simulations/outputs/longsims/longsims_n20_a035_lowvar.csv')
  )

### Get summaries of rescue for each trial

# Get data on size-rescue (three successive generations exceeding original size)
rescue.summary = all.data %>%
  # Get generations *after initial* where population size exceeds original
  filter(n > n.pop0) %>%
  # For each trial, get the generation in which this first occurs
  group_by(trial, n.pop0, low.var, alpha) %>%
  summarise(rescue.time = which.diff(gen, 3)) %>%
  # Get only trials which are rescued
  # (because which.diff() returns 0 if there is no such event)
  filter(rescue.time > 0) %>%
  ungroup()

# Get overall probability of rescue
p.rescue = merge(
  # Get initial genotype for each trial
  x = all.data %>% 
    distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
    select(trial, n.pop0, low.var, alpha, gbar),
  # Merge this with the rescue summary from above
  y = rescue.summ,
  by = c('trial', 'n.pop0', 'low.var', 'alpha'),
  all.x = TRUE, all.y = TRUE
) %>%
  # `rescued` column is response - NAs for trials with no designed rescue time
  # (n.b. these rows were removed above)
  mutate(rescued = !is.na(rescue.time))

# Get data on time of fitness-based rescue (three successive generations with intrinsic fitness > 1)
wescue.summary = all.data %>%
  # Define "extrinsic fitness", expected reproductive output
  # (intrinsic fitness times density dependent term)
  mutate(rbar = wbar * exp(-alpha * n)) %>%
  # Get only cases (after first generation) where fitness exceeds 1
  filter(rbar > 1) %>%
  # Get the time in which this occurs
  group_by(trial, n.pop0, low.var, alpha) %>%
  summarise(rescue.time = which.diff(gen, 3)) %>%
  # Remove cases where rescue has not occurred
  filter(rescue.time > 0)

# Get overall probability of rescue
p.wescue = merge(
  # Get initial genotype for each trial
  x = all.data %>% 
    distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
    select(trial, n.pop0, low.var, alpha, gbar),
  # Merge this with the rescue summary from above
  y = wescue.summary,
  by = c('trial', 'n.pop0', 'low.var', 'alpha'),
  all.x = TRUE, all.y = TRUE
) %>%
  # `rescued` column is response - NAs for trials with no designed rescue time
  # (n.b. these rows were removed above)
  mutate(rescued = !is.na(rescue.time))

### Run statistical models

# Formatting
p.rescue = p.rescue %>%
  mutate(size = factor(ifelse(n.pop0 > 50, 'Large', 'Small')),
         gdiv = factor(ifelse(low.var, 'Low', 'High')),
         ddep = factor(ifelse(alpha == 0, 'Ind', 'Neg')))

rescue.mod.full =  stan_glm(
  rescued ~ size * gdiv * ddep * gbar,
  data = p.rescue, 
  seed = 89210,
  cores = 4,
  prior = normal(0, 5),
  family = 'binomial'
)

summary(rescue.mod.full)

# Posterior check
bayesplot::pp_check(rescue.mod.full)
# ccc = posterior_predict(rescue.mod.full)
# cc2 = apply(ccc, 2, mean)
# hist(cc2)

rescue.post.df = as.data.frame(rescue.mod.full)

# NDD effects in high diversity, large populations
with(rescue.post.df, mean(ddepNeg))

# NDD effecrts in high diversity, small populations
with(rescue.post.df, mean(ddepNeg + `sizeSmall:ddepNeg`))

# NDD effects in low diversity, large populations
with(rescue.post.df, mean(ddepNeg + `gdivLow:ddepNeg`))

# NDD effects in low diversity, small populations
with(rescue.post.df, mean(ddepNeg + `sizeSmall:ddepNeg` + `gdivLow:ddepNeg` + `sizeSmall:gdivLow:ddepNeg`))

## Effect sizes for density dependence
ndd.rescue = with(
  rescue.post.df,
  cbind(
    hidiv.large = ddepNeg,
    hidiv.small = ddepNeg + `sizeSmall:ddepNeg`,
    lodiv.large = ddepNeg + `gdivLow:ddepNeg`,
    lodiv.small = ddepNeg + `sizeSmall:ddepNeg` + `gdivLow:ddepNeg` + `sizeSmall:gdivLow:ddepNeg`
  )
)

apply(ndd.rescue, 2, mean)
apply(ndd.rescue, 2, function(x) c(quantile(x, 0.025), quantile(x, 0.975)))

## Effect sizes for size (being small)
size.rescue = with(
  rescue.post.df,
  cbind(
    hidiv.did = sizeSmall,
    hidiv.ndd = sizeSmall + `sizeSmall:ddepNeg`,
    lodiv.did = sizeSmall + `sizeSmall:gdivLow`,
    lodiv.ndd = sizeSmall + `sizeSmall:ddepNeg` + `sizeSmall:gdivLow` + `sizeSmall:gdivLow:ddepNeg`
  )
)

apply(size.rescue, 2, mean)
apply(size.rescue, 2, function(x) c(quantile(x, 0.025), quantile(x, 0.975)))

## Effect sizes for low diversity 
lodiv.rescue = with(
  rescue.post.df,
  cbind(
    large.did = gdivLow,
    large.ndd = gdivLow + `gdivLow:ddepNeg`,
    small.did = gdivLow + `sizeSmall:gdivLow`,
    small.ndd = gdivLow + `gdivLow:ddepNeg` + `sizeSmall:gdivLow` + `sizeSmall:gdivLow:ddepNeg`
  )
)

apply(lodiv.rescue, 2, mean)
apply(lodiv.rescue, 2, function(x) c(quantile(x, 0.025), quantile(x, 0.975)))

## Effect size of increasing maladaptation
gbar.slopes = with(
  -1*rescue.post.df,
  cbind(
    large.hidiv.did = gbar,
    large.hidiv.ndd = gbar + `ddepNeg:gbar`,
    large.lodiv.did = gbar + `gdivLow:gbar`,
    large.lodiv.ndd = gbar + `ddepNeg:gbar` + `gdivLow:gbar` + `gdivLow:ddepNeg:gbar`,
    small.hidiv.did = gbar + `sizeSmall:gbar`,
    small.hidiv.ndd = gbar + `sizeSmall:gbar` + `ddepNeg:gbar` + `sizeSmall:ddepNeg:gbar`,
    small.lodiv.did = gbar + `sizeSmall:gbar` + `gdivLow:gbar` + `sizeSmall:gdivLow:gbar`,
    small.lodiv.ndd = gbar + `sizeSmall:gbar` + `gdivLow:gbar` + `ddepNeg:gbar` + 
                      `sizeSmall:gdivLow:gbar` + `sizeSmall:ddepNeg:gbar` + `gdivLow:ddepNeg:gbar` +
                      `sizeSmall:gdivLow:ddepNeg:gbar`
  )
)

apply(gbar.slopes, 2, mean)
apply(gbar.slopes, 2, function(x) c(quantile(x, 0.025), quantile(x, 0.975)))

### Wescue (fitness rescue) models

# Formatting
p.wescue = p.wescue %>%
  mutate(size = factor(ifelse(n.pop0 > 50, 'Large', 'Small')),
         gdiv = factor(ifelse(low.var, 'Low', 'High')),
         ddep = factor(ifelse(alpha == 0, 'Ind', 'Neg')))

wescue.mod.full =  stan_glm(
  rescued ~ size * gdiv * ddep * gbar,
  data = p.wescue, 
  seed = 422700,
  cores = 4,
  prior = normal(0, 5),
  family = 'binomial'
)

summary(wescue.mod.full)

bayesplot::pp_check(wescue.mod.full)

# (next: effect sizes)