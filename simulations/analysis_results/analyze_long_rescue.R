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
  y = rescue.summary,
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

apply(ndd.rescue, 2, mean) %>% round(2)
apply(ndd.rescue, 2, function(x) c(quantile(x, 0.025), quantile(x, 0.975))) %>% round(2)

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

apply(size.rescue, 2, mean) %>% round(2)
apply(size.rescue, 2, function(x) c(quantile(x, 0.025), quantile(x, 0.975))) %>% round(2)

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

apply(lodiv.rescue, 2, mean) %>% round(2)
apply(lodiv.rescue, 2, function(x) c(quantile(x, 0.025), quantile(x, 0.975))) %>% round(2)

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

apply(gbar.slopes, 2, mean) %>% round(2)
apply(gbar.slopes, 2, function(x) c(quantile(x, 0.025), quantile(x, 0.975))) %>% round(2)

### Wescue (fitness rescue) models

# Formatting
p.wescue = p.wescue %>%
  mutate(size = factor(ifelse(n.pop0 > 50, 'Large', 'Small')),
         gdiv = factor(ifelse(low.var, 'Low', 'High')),
         ddep = factor(ifelse(alpha == 0, 'Ind', 'Neg')))

# Run model
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

wescue.post.df = as.data.frame(wescue.mod.full)

# Effects of density dependence
ndd.wescue = with(
  wescue.post.df,
  cbind(
    hidiv.large = ddepNeg,
    hidiv.small = ddepNeg + `sizeSmall:ddepNeg`,
    lodiv.large = ddepNeg + `gdivLow:ddepNeg`,
    lodiv.small = ddepNeg + `sizeSmall:ddepNeg` + `gdivLow:ddepNeg` + `sizeSmall:gdivLow:ddepNeg`
  )
)

apply(ndd.wescue, 2, mean) %>% round(2)
apply(ndd.wescue, 2, function(x) c(quantile(x, 0.025), quantile(x, 0.975))) %>% round(2)

# Effects of size
size.wescue = with(
  wescue.post.df,
  cbind(
    hidiv.did = sizeSmall,
    hidiv.ndd = sizeSmall + `sizeSmall:ddepNeg`,
    lodiv.did = sizeSmall + `sizeSmall:gdivLow`,
    lodiv.ndd = sizeSmall + `sizeSmall:ddepNeg` + `sizeSmall:gdivLow` + `sizeSmall:gdivLow:ddepNeg`
  )
)

apply(size.wescue, 2, mean) %>% round(2)
apply(size.wescue, 2, function(x) c(quantile(x, 0.025), quantile(x, 0.975))) %>% round(2)

# Effects of low diversity
lodiv.wescue = with(
  wescue.post.df,
  cbind(
    large.did = gdivLow,
    large.ndd = gdivLow + `gdivLow:ddepNeg`,
    small.did = gdivLow + `sizeSmall:gdivLow`,
    small.ndd = gdivLow + `gdivLow:ddepNeg` + `sizeSmall:gdivLow` + `sizeSmall:gdivLow:ddepNeg`
  )
)

apply(lodiv.wescue, 2, mean) %>% round(2)
apply(lodiv.wescue, 2, function(x) c(quantile(x, 0.025), quantile(x, 0.975))) %>% round(2)

# Effects of genotype (increased adaptation)
gbar.slopes.w = with(
  -1*wescue.post.df,
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

apply(gbar.slopes.w, 2, mean) %>% round(2)
apply(gbar.slopes.w, 2, function(x) c(quantile(x, 0.025), quantile(x, 0.975))) %>% round(2)

### Plots

# Posterior draws from rescue model
geno.res.fulls = expand.grid(gbar = (-5:5)/10,
                             ddep = c("Ind", "Neg"),
                             size = c("Large", "Small"),
                             gdiv = c("Low", "High"))

res.epreds = posterior_epred(rescue.mod.full, 
                             newdata = geno.res.fulls,
                             seed = 8400, draws = 200) %>%
  as.data.frame() %>%
  t() %>%
  cbind(geno.res.fulls, .) %>%
  gather(key = draw, value = estimate, -c(gbar, ddep, size, gdiv))

# Generate plot of Pr(rescue) ~ genotype and other factors
# for rescue based on population size
res.geno.plot = res.epreds %>%
  mutate(size = paste("Initially", tolower(size)),
         gdiv = paste(gdiv, "diversity")) %>%
  ggplot(aes(x = 2.75 - gbar, y = estimate)) +
  geom_line(
    aes(
      group = interaction(ddep, draw),
      colour = factor(ddep)
    ),
    size = 0.1
  ) +
  labs(x = '', y = 'Probability of rescue\n(exceeding initial size)') +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap( ~ paste(size, gdiv, sep = ', '), ncol = 4) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title  = element_text(size = 8),
        strip.text  = element_text(size = 8),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        plot.margin = margin(b = 0, r = 5, l = 20, unit = 'pt'))

res.geno.plot

# Generate posterior draws from rescue model (based on fitness)
geno.wes.fulls = expand.grid(gbar = (-5:5)/10,
                             ddep = c("Ind", "Neg"),
                             size = c("Large", "Small"),
                             gdiv = c("Low", "High"))

wes.epreds = posterior_epred(wescue.mod.full, 
                             newdata = geno.wes.fulls,
                             seed = 8400, draws = 200) %>%
  as.data.frame() %>%
  t() %>%
  cbind(geno.wes.fulls, .) %>%
  gather(key = draw, value = estimate, -c(gbar, ddep, size, gdiv))


wes.geno.plot = ggplot(wes.epreds, aes(x = 2.75 - gbar, y = estimate)) +
  geom_line(
    aes(
      group = interaction(ddep, draw),
      colour = factor(ddep)
    ),
    size = 0.1
  ) +
  labs(x = 'Initial maladaptation', y = 'Probability of rescue\n(mean fitness exceeding 1)') +
  scale_color_manual(values = c('black', 'purple')) +
  facet_wrap( ~ paste(size, gdiv, sep = ', '), ncol = 4) +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 8, angle = 45),
        axis.text.y = element_text(size = 8),
        axis.title  = element_text(size = 8),
        strip.text  = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        plot.margin = margin(b = 0, r = 5, l = 20, unit = 'pt'))

wes.geno.plot

geno.plot.grid = plot_grid(res.geno.plot, 
                           wes.geno.plot,
                           labels = c("A)", "B)"),
                           nrow = 2)

save_plot(geno.plot.grid, base_width = 8, base_height = 5,
          filename = 'simulations/analysis_results/figure_drafts/draft_figs/fig_supp_p_rescue.pdf')

### Next: instant rescue

rescue.summary %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(inst.rescue = sum(rescue.time < 3),
            total.rescue = n(),
            propn.rescue = mean(rescue.time < 3))

wescue.summary %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(inst.rescue = sum(rescue.time < 3),
            total.rescue = n(),
            propn.rescue = mean(rescue.time < 3))
