# Script assessing probability of rescue in long (50-generation) simulations

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
  if(length(x) >= lag) { 
    if (any(diff(x, lag = (lag - 1)) == lag - 1)) { x[which.max(diff(x, lag = (lag - 1)) == (lag-1))] 
    } else { NA }
  } else { NA }
}

# Load in simulation data
all.data = 
  rbind(
    read.csv('simulations/outputs/longsims_n100_a000_hivar.csv'),
    read.csv('simulations/outputs/longsims_n100_a000_lowvar.csv'),
    read.csv('simulations/outputs/longsims_n20_a000_hivar.csv'),
    read.csv('simulations/outputs/longsims_n20_a000_lowvar.csv'),
    read.csv('simulations/outputs/longsims_n100_a035_hivar.csv'),
    read.csv('simulations/outputs/longsims_n100_a035_lowvar.csv'),
    read.csv('simulations/outputs/longsims_n20_a035_hivar.csv'),
    read.csv('simulations/outputs/longsims_n20_a035_lowvar.csv')
  ) %>%
  # Re-index generation starting at 1
  mutate(gen = gen - 1)

### Get summaries of rescue for each trial (for both types of rescue)

# Fitness-based rescue ("w"escue because fitness... w)

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
  filter(!is.na(rescue.time))

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
  # (n.b. these rows were removed above when defining `wescue.summary`)
  mutate(rescued = !is.na(rescue.time))

head(p.wescue)

### Next: size-based rescue

# Get data on size-rescue (three successive generations exceeding original size)
rescue.summary = all.data %>%
  # Get generations *after initial* where population size exceeds original
  filter(n > n.pop0) %>%
  # For each trial, get the generation in which this first occurs
  group_by(trial, n.pop0, low.var, alpha) %>%
  summarise(rescue.time = which.diff(gen, 3)) %>%
  # Get only trials which are rescued
  # (because which.diff() returns 0 if there is no such event)
  filter(!is.na(rescue.time)) %>%
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

##### Information for Tables SH1 and SH2s

# Probabilities of rescue (overall) and proportion of rescues that are
# instantaneous

p.wescue %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(p.rescue = round(mean(rescued), 2)) %>%
  arrange(desc(n.pop0, low.var, alpha))

wescue.summary %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(inst.rescue = sum(rescue.time < 2),
            total.rescue = n(),
            propn.rescue = round(mean(rescue.time < 2), 2)) %>%
  arrange(desc(n.pop0, low.var, alpha))

p.rescue %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(p.rescue = round(mean(rescued), 2)) %>%
  arrange(desc(n.pop0, low.var, alpha))

rescue.summary %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(inst.rescue = sum(rescue.time < 2),
            total.rescue = n(),
            propn.rescue = round(mean(rescue.time < 2), 2)) %>%
  arrange(desc(n.pop0, low.var, alpha))

### Times to rescue 

# Get mean times to fitness-rescue
wescue.time.means = wescue.summary %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(mean.t = mean(rescue.time)) %>%
  ungroup() %>%
  mutate(size = paste("Initially", ifelse(n.pop0 > 50, 'large', 'small')),
         gdiv = paste(ifelse(low.var, 'Low', 'high'), "diversity"))

wescue.time.means %>%
  arrange(desc(n.pop0, low.var, alpha)) %>%
  mutate(mean.t = round(mean.t, 1))

# Get mean times to size-rescue
rescue.time.means = rescue.summary %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(mean.t = mean(rescue.time)) %>%
  ungroup() %>%
  mutate(size = paste("Initially", ifelse(n.pop0 > 50, 'large', 'small')),
         gdiv = paste(ifelse(low.var, 'Low', 'high'), "diversity"))

rescue.time.means %>%
  arrange(desc(n.pop0, low.var, alpha)) %>%
  mutate(mean.t = round(mean.t, 1))

# Comparison of effects of NDD
wescue.time.means %>%
  spread(key = alpha, value = mean.t)

rescue.time.means %>%
  spread(key = alpha, value = mean.t)

### Extinctions after rescue

# Get which populations went extinct
n.extinct = all.data %>%
  # Get final generation's record for each population
  group_by(trial, n.pop0, low.var, alpha)  %>%
  slice_max(gen) %>%
  # Add "extinct" denotation
  mutate(extinct = gen < 50 & n < 5000) %>%
  select(trial, gen, n, n.pop0, low.var, alpha, extinct)

# Merge extinct populations with the fitness-based rescue list
# (all.y = FALSE means only rescued populations will be in this data frame)
merge(wescue.summary, n.extinct, all.y = FALSE) %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(p.extinct = mean(extinct)) %>%
  arrange(desc(n.pop0, low.var, alpha)) %>%
  mutate(p.extinct = round(p.extinct, 2))

merge(rescue.summary, n.extinct, all.y = FALSE) %>%
  group_by(n.pop0, low.var, alpha) %>%
  summarise(p.extinct = mean(extinct)) %>%
  arrange(desc(n.pop0, low.var, alpha)) %>%
  mutate(p.extinct = round(p.extinct, 2))

### Run statistical models in Rstanarm

### Wescue (fitness rescue) models

# Formatting (strings for readability)
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
# Rhat < 0, some small n_eff but should be sufficiently large for estimates

# Posterior check
bayesplot::pp_check(wescue.mod.full)

# Convert to data frame (to get posterior estimates)
wescue.post.df = as.data.frame(wescue.mod.full)

# Looking at effect sizes

# Effects of density dependence
ndd.wescue = with(
  wescue.post.df,
  cbind(
    hidiv.large = ddepNeg,
    lodiv.large = ddepNeg + `gdivLow:ddepNeg`,
    hidiv.small = ddepNeg + `sizeSmall:ddepNeg`,
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

### Assessment of size-based rescue

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

# Get posteriors from data frame
rescue.post.df = as.data.frame(rescue.mod.full)

## Effect sizes for density dependence
ndd.rescue = with(
  rescue.post.df,
  cbind(
    hidiv.large = ddepNeg,
    lodiv.large = ddepNeg + `gdivLow:ddepNeg`,
    hidiv.small = ddepNeg + `sizeSmall:ddepNeg`,
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

### Figure F1 (supplements) - probabilities of rescue as function of genotype

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
  # mutate(size = paste("Initially", tolower(size)),
  #        gdiv = paste(gdiv, "diversity")) %>%
  ggplot(aes(x = 2.8 - gbar, y = estimate)) +
  geom_line(
    aes(
      group = interaction(ddep, draw),
      colour = factor(ddep)
    ),
    size = 0.1
  ) +
  labs(x = 'Initial maladaptation', 
       y = 'Probability of size-rescue') +
  scale_color_manual(values = c('black', 'purple')) +
  scale_x_continuous(breaks = 2 * (12:16) / 10) +
  facet_wrap( ~ paste(size, gdiv, sep = ', '), ncol = 4) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title  = element_text(size = 8),
        strip.text  = element_blank(),
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


wes.geno.plot = wes.epreds %>%
  mutate(size = paste("Initially", tolower(size)),
         gdiv = paste(gdiv, "diversity")) %>%
  ggplot(aes(x = 2.8 - gbar, y = estimate)) +
  geom_line(
    aes(
      group = interaction(ddep, draw),
      colour = factor(ddep)
    ),
    size = 0.1
  ) +
  labs(x = '', y = 'Probability of fitness-rescue') +
  scale_color_manual(values = c('black', 'purple')) +
  scale_x_continuous(breaks = 2 * (12:16) / 10) +
  facet_wrap( ~ paste(size, gdiv, sep = ', '), ncol = 4) +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 8, angle = 45),
        axis.text.y = element_text(size = 8),
        axis.title  = element_text(size = 8),
        strip.text  = element_text(size = 8),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        plot.margin = margin(b = 0, r = 5, l = 20, unit = 'pt'))

wes.geno.plot

geno.plot.grid = plot_grid(wes.geno.plot, 
                           res.geno.plot,
                           labels = c("A)", "B)"),
                           nrow = 2)

save_plot(geno.plot.grid, base_width = 8, base_height = 5,
          filename = 'analysis_results/figures/fig_h1_p_rescue.png')

### Time until rescue figure (supporting info)

# Time until fitness-rescue histogram
wescue.time.hist = wescue.summary %>% 
  ungroup() %>%
  mutate(rescue.time = rescue.time - min(rescue.time) + 1) %>%
  mutate(size = paste("Initially", ifelse(n.pop0 > 50, 'large', 'small')),
         gdiv = paste(ifelse(low.var, 'Low', 'high'), "diversity")) %>%
  ggplot(aes(x = rescue.time)) +
  geom_histogram(
    aes(
      fill = factor(alpha)
    ),
    position = 'identity',
    binwidth = 1, 
    alpha = 0.5) +
  geom_point(
    data = wescue.time.means,
    aes(
      x = mean.t, 
      y = 240,
      colour = factor(alpha) 
    ),
    size = 4, shape = 18
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  facet_wrap(~ paste(size, gdiv, sep = ', '), nrow = 1) +
  labs(x = '', y = 'Frequency\n(fitness rescue)') +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )

# Time until size-rescue histogram
rescue.time.hist = rescue.summary %>% 
  ungroup() %>%
  mutate(rescue.time = rescue.time - min(rescue.time) + 1) %>%
  mutate(size = paste("Initially", ifelse(n.pop0 > 50, 'large', 'small')),
         gdiv = paste(ifelse(low.var, 'Low', 'high'), "diversity")) %>%
  ggplot(aes(x = rescue.time)) +
  geom_histogram(
    aes(
      fill = factor(alpha)
    ),
    position = 'identity',
    binwidth = 2, 
    alpha = 0.5) +
  geom_point(
    data = rescue.time.means,
    aes(
      x = mean.t, 
      y = 220,
      colour = factor(alpha) 
    ),
    size = 4, shape = 18
  ) +
  scale_color_manual(values = c('black', 'purple')) +
  scale_fill_manual(values = c('black', 'purple')) +
  labs(x = 'Generation', y = 'Frequency\n(size rescue)') +
  facet_wrap(~ paste(size, gdiv, sep = ', '), nrow = 1) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = 'none',
    panel.background = element_blank()
  )

# Combine histograms
time.hists = plot_grid(
  wescue.time.hist, rescue.time.hist,
  labels = c("A)", "B)"),
  nrow = 2
)

time.hists

ggsave('analysis_results/figures/fig_h2_rescue_times.png',
       width = 8, height = 5)

### Session info
# (as of October 24 2022)

# R version 4.1.0 (2021-05-18)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 11.2.3
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] tidyr_1.1.3     dplyr_1.0.7     cowplot_1.1.1   ggplot2_3.3.4   rstanarm_2.21.1 Rcpp_1.0.7     
# 
# loaded via a namespace (and not attached):
#   [1] splines_4.1.0        jsonlite_1.7.2       gtools_3.9.2         RcppParallel_5.1.4   StanHeaders_2.21.0-7
# [6] threejs_0.3.3        shiny_1.6.0          stats4_4.1.0         pillar_1.6.1         lattice_0.20-44     
# [11] glue_1.4.2           digest_0.6.27        promises_1.2.0.1     minqa_1.2.4          colorspace_2.0-2    
# [16] Matrix_1.3-3         htmltools_0.5.1.1    httpuv_1.6.2         plyr_1.8.6           dygraphs_1.1.1.6    
# [21] pkgconfig_2.0.3      rstan_2.21.2         purrr_0.3.4          xtable_1.8-4         scales_1.1.1        
# [26] processx_3.5.2       later_1.3.0          lme4_1.1-27.1        tibble_3.1.2         bayesplot_1.8.1     
# [31] generics_0.1.0       ellipsis_0.3.2       DT_0.18              withr_2.4.2          shinyjs_2.0.0       
# [36] cli_2.5.0            survival_3.2-11      magrittr_2.0.1       crayon_1.4.1         mime_0.11           
# [41] ps_1.6.0             fansi_0.5.0          nlme_3.1-152         MASS_7.3-54          xts_0.12.1          
# [46] pkgbuild_1.2.0       colourpicker_1.1.0   rsconnect_0.8.24     tools_4.1.0          loo_2.4.1           
# [51] prettyunits_1.1.1    lifecycle_1.0.0      matrixStats_0.60.1   stringr_1.4.0        V8_3.4.2            
# [56] munsell_0.5.0        callr_3.7.0          compiler_4.1.0       rlang_0.4.11         nloptr_1.2.2.2      
# [61] grid_4.1.0           ggridges_0.5.3       rstudioapi_0.13      htmlwidgets_1.5.3    crosstalk_1.1.1     
# [66] igraph_1.2.6         miniUI_0.1.1.1       base64enc_0.1-3      boot_1.3-28          gtable_0.3.0        
# [71] codetools_0.2-18     inline_0.3.19        DBI_1.1.1            curl_4.3.2           markdown_1.1        
# [76] reshape2_1.4.4       R6_2.5.0             gridExtra_2.3        rstantools_2.1.1     zoo_1.8-9           
# [81] fastmap_1.1.0        utf8_1.2.1           shinystan_2.5.0      shinythemes_1.2.0    stringi_1.7.3       
# [86] parallel_4.1.0       vctrs_0.3.8          tidyselect_1.1.1    
