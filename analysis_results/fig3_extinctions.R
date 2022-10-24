# Script containing all results pertaining to extinctions in the long (50
# generation) simulations.
# Contains code for Figure 2
# SN - 29 Mar 2021
# re-run for archiving 4 Oct 2022

# Load packages
library(rstanarm)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)

### Setup and declaring variables

# Clear namespace
rm(list = ls())

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
  # Re-index time to start at generation 0 instead of 1
  mutate(gen = gen - 1)

# Add in extinction designation
all.data = all.data %>%
  # Get final generation's record for each population
  group_by(trial, n.pop0, low.var, alpha)  %>%
  slice_max(gen) %>%
  # Add "extinct" denotation
  group_by(n.pop0, low.var, alpha, trial) %>%
  summarise(extinct = gen < 50 & n < 5000) %>%
  # Combine this "extinct" column with main data frame
  merge(x = all.data, y = ., by = c('trial', 'n.pop0', 'low.var', 'alpha'))

# Get number of extinctions occurring in each generation
all.extinctions.long = all.data %>%
  select(gen, trial, n.pop0, low.var, alpha, extinct) %>%
  group_by(trial, n.pop0, low.var, alpha) %>%
  mutate(ext.this.gen = (gen %in% max(gen)) & extinct) %>%
  group_by(n.pop0, low.var, alpha, gen) %>%
  summarise(n.extinctions = sum(ext.this.gen)) %>%
  ungroup() %>%
  arrange(n.pop0, low.var, alpha, gen) %>%
  # Reformating levels into plot-readible labels
  mutate(romnum = factor(ifelse(n.pop0 %in% 100,
                                ifelse(low.var, '(i)', '(ii)'),
                                ifelse(low.var, '(iii)', '(iv)'))),
         n.pop0 = factor(ifelse(n.pop0 %in% 100, "Large", "Small")),
         alpha = factor(ifelse(alpha > 0, "Density dependent", "Density independent")),
         low.var = factor(ifelse(low.var, "Low diversity", "High diversity")))

##### Bayesian model for Figure 3D and effect sizes

# Get data for extinction models:
# one row (entry) per simulation - using the first generation here to capture
# initial genotype.
ext.data = all.data %>% 
  distinct(trial, n.pop0, low.var, alpha, .keep_all = TRUE) %>%
  mutate(
    low.var = factor(ifelse(low.var, 'low.div', 'hi.div')),
    n.pop0  = factor(ifelse(n.pop0 > 50, 'large', 'small')),
    alpha   = factor(ifelse(alpha > 0, 'ndd', 'di'))
  )

# Fit model
ext.mod.full =  stan_glm(
  extinct ~ n.pop0 * low.var * alpha * gbar,
  data = ext.data, 
  seed = 992,
  cores = 4,
  prior = normal(0, 5),
  family = 'binomial'
)

summary(ext.mod.full)

# Model checks
bayesplot::pp_check(ext.mod.full)
plot(ext.mod.full$fitted.values)

# Posterior estimates for getting effect sizes
posts = as.data.frame(ext.mod.full)

### Effect sizes for Table 2

# Effect sizes for density dependence for all combos of other treatments

ndd.posts = with(
  posts,
  cbind(
    large.hivar = alphandd,
    large.lovar = alphandd + `low.varlow.div:alphandd`,
    small.hivar = alphandd + `n.pop0small:alphandd`,
    small.lovar = alphandd + `low.varlow.div:alphandd` + `n.pop0small:alphandd` + `n.pop0small:low.varlow.div:alphandd`
  )
)

# Mean and 95% CI
apply(ndd.posts, 2, mean) %>% round(2)
apply(ndd.posts, 2, function(x) quantile(x, c(0.025, 0.975))) %>% round(2)

# Effect sizes for size (small) for all combos of other treatments

size.posts = with(
  posts,
  cbind(
    hivar.dindp = n.pop0small,
    hivar.ddepn = n.pop0small + `n.pop0small:alphandd`,
    lovar.dindp = n.pop0small + `n.pop0small:low.varlow.div`,
    lovar.ddepn = n.pop0small + `n.pop0small:alphandd` + `n.pop0small:low.varlow.div` + `n.pop0small:low.varlow.div:alphandd`
  )
)

apply(size.posts, 2, mean) %>% round(2)
apply(size.posts, 2, function(x) quantile(x, c(0.025, 0.975))) %>% round(2)

# Effect sizes for low genetic diversity for all combos of other treatments

gdiv.posts = with(
  posts,
  cbind(
    large.dindp = low.varlow.div,
    large.ddepn = low.varlow.div + `low.varlow.div:alphandd`,
    small.dindp = low.varlow.div + `n.pop0small:low.varlow.div`,
    small.ddepn = low.varlow.div + `low.varlow.div:alphandd` + `n.pop0small:low.varlow.div` + `n.pop0small:low.varlow.div:alphandd`
  )
)

apply(gdiv.posts, 2, mean) %>% round(2)
apply(gdiv.posts, 2, function(x) quantile(x, c(0.025, 0.975))) %>% round(2)

# Effect s izes for increasing maladaptation

gbar.slopes = with(
  -1*posts,
  cbind(
    large.hivar.did = gbar,
    small.hivar.did = gbar + `n.pop0small:gbar`,
    large.lovar.did = gbar + `low.varlow.div:gbar`,
    large.hivar.ndd = gbar + `alphandd:gbar`,
    small.lovar.did = gbar + `n.pop0small:gbar` + `low.varlow.div:gbar` + `n.pop0small:low.varlow.div:gbar`,
    small.hivar.ndd = gbar + `n.pop0small:gbar` + `alphandd:gbar` + `n.pop0small:alphandd:gbar`,
    large.lovar.ndd = gbar + `low.varlow.div:gbar` + `alphandd:gbar` + `low.varlow.div:alphandd:gbar`,
    small.lovar.ndd = gbar + `n.pop0small:gbar` + `low.varlow.div:gbar` + `alphandd:gbar` +
      `n.pop0small:low.varlow.div:gbar` + `n.pop0small:alphandd:gbar` + `low.varlow.div:alphandd:gbar` + 
      `n.pop0small:low.varlow.div:alphandd:gbar`
  )
)

apply(gbar.slopes, 2, mean) %>% round(2)
apply(gbar.slopes, 2, function(x) quantile(x, c(0.025, 0.975))) %>% round(2)

# Arranged to highlight NDD effects
apply(gbar.slopes, 2, mean)[c(1,4,3,7,2,6,5,8)] %>% round(2)
apply(gbar.slopes, 2, function(x) quantile(x, c(0.025, 0.975)))[,c(1,4,3,7,2,6,5,8)] %>% round(2)

### Other summary statistics

# Increased extinction risk for each parameter combo
ext.data %>% 
  group_by(n.pop0, low.var, alpha) %>%
  summarise(pext = mean(extinct)) %>%
  spread(key = alpha, value = pext) %>%
  mutate(incr = ndd / di)

### Try a figure with opposite direction of paneling

size.plot.row = all.data %>%
  filter(trial < 26) %>%
  mutate(romnum = factor(ifelse(n.pop0 %in% 100,
                                ifelse(!low.var, '(i)', '(ii)'),
                                ifelse(!low.var, '(iii)', '(iv)'))),
         n.pop0 = factor(ifelse(n.pop0 %in% 100, "Large", "Small")),
         alpha = factor(ifelse(alpha > 0, "Density dependent", "Density independent")),
         low.var = factor(ifelse(!low.var, "High diversity", "Low diversity"))) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      x = gen, 
      y = n, 
      colour = factor(alpha), 
      linetype = extinct, 
      group = interaction(alpha, trial)
    ),
    size = 0.6) +
  geom_line(
    aes(
      x = gen, 
      y = n, 
      colour = factor(alpha), 
      linetype = extinct, 
      group = interaction(alpha, trial),
      alpha = factor(alpha)
    ),
    size = 0.5) +
  labs(x = 'Generation', y = '') +
  scale_y_log10() +
  scale_alpha_manual(values = c(0, 1)) +
  scale_fill_manual(values = c('purple', 'black')) +
  scale_color_manual(values = c('purple', 'black')) +
  facet_wrap( ~ paste0(romnum, ' ', n.pop0, ', ', low.var), ncol = 4) +
  theme(legend.position = 'none',
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        strip.background = element_rect(colour = 'black'),
        strip.text = element_text(size = 11))

inst.plot.row = all.extinctions.long %>%
  group_by(n.pop0, low.var, alpha) %>%
  mutate(n.extant = 1000 - c(0, cumsum(n.extinctions)[-50]),
         p.instant.extinct = n.extinctions / n.extant) %>%
  ggplot(aes(x = gen)) +
  geom_ribbon(
    aes(
      ymin = p.instant.extinct - 2 * sqrt(p.instant.extinct*(1-p.instant.extinct) / n.extant),
      ymax = p.instant.extinct + 2 * sqrt(p.instant.extinct*(1-p.instant.extinct) / n.extant),
      fill = factor(alpha)
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      y = p.instant.extinct,
      group = factor(alpha),
      colour = factor(alpha)
    ),
    size = 1
  ) +
  labs(x = 'Generation', y = '') +
  scale_fill_manual(values = c('purple', 'black')) +
  scale_color_manual(values = c('purple', 'black')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        plot.margin = margin(b = 0, r = 5, l = 20, unit = 'pt'))

# Cumulative extinction plot (b)
cuml.plot.row = all.extinctions.long %>%
  group_by(n.pop0, low.var, alpha) %>%
  mutate(p.extinct = cumsum(n.extinctions / 1000)) %>%
  ggplot(aes(x = gen)) +
  geom_ribbon(
    aes(
      xmin = 0, xmax = 50,
      ymin = 0, ymax = p.extinct,
      group = interaction(n.pop0, low.var, alpha),
      fill = factor(alpha)
    ),
    alpha = 0.5
  ) +
  geom_line(
    aes(
      y = p.extinct,
      group = factor(alpha),
      colour = factor(alpha)
    )
  ) +
  labs(x = 'Generation', y = '') +
  scale_fill_manual(values = c('purple', 'black')) +
  scale_colour_manual(values = c('purple', 'black')) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  theme(legend.position = 'none',        
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        plot.margin = margin(b = 0, r = 5, l = 20, unit = 'pt'))

### Get estimates for Fig. 3D (posterior estimates of extinction probability
#   according to GLM)

# Get "genotypes" for which to get 
geno.ext.fulls = expand.grid(
  gbar = (-5:5)/10,
  alpha = c("di", "ndd"),
  n.pop0 = c("large", "small"),
  low.var = c("low.div", "hi.div")
) %>%
  filter(n.pop0 %in% 'small' | abs(gbar) < 0.25)

epreds = posterior_epred(
  ext.mod.full, 
  newdata = geno.ext.fulls,
  seed = 196200, 
  draws = 200
) %>%
  as.data.frame() %>%
  t() %>%
  cbind(geno.ext.fulls, .)

epreds = epreds %>%
  gather(key = draw, value = estimate, -c(gbar, alpha, n.pop0, low.var))

geno.full.row = ggplot(epreds, aes(x = 2.8 - gbar, y = estimate)) +
  geom_line(
    aes(
      group = interaction(alpha, draw),
      colour = factor(alpha)
    ),
    size = 0.1
  ) +
  geom_rug(
    data = ext.data,
    inherit.aes = FALSE,
    aes(
      x = 2.8 - gbar
    ),
    sides = "b",
    size = 0.1
  ) +
  labs(x = 'Initial maladaptation', y = '') +
  scale_color_manual(values = c('black', 'purple')) +
  scale_x_continuous(breaks = 2 * (12:16) / 10) +
  facet_wrap( ~ paste(n.pop0, low.var, sep = ', '), ncol = 4) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = 'gray88'),
        panel.background = element_rect(fill = 'white'),
        plot.margin = margin(b = 0, r = 5, l = 20, unit = 'pt'))

geno.full.row

data.plots.row = plot_grid(size.plot.row, inst.plot.row, 
                           cuml.plot.row, geno.full.row, 
                           labels = c('(A)', '(B)', '(C)', '(D)'),
                           label_x = c(-0.01, -0.01, -0.01, -0.01),
                           label_size = 12,
                           nrow = 4)

data.plots.row

# Get legend for figure
extinct.legend = get_legend(
  inst.plot.row + 
    guides(fill = guide_legend(''), colour = guide_legend('')) +
    theme(legend.position = 'bottom',
          legend.box.margin = margin(1, 0, 0, 0))
)

# Export plot
plot_grid(data.plots.row, extinct.legend, ncol = 1, rel_heights = c(1, .1)) %>%
  save_plot(filename = 'analysis_results/figures/fig_ext_fig.png',
            base_width = 8, base_height = 8)

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
