# Script for simulating population dynamics of an adapting bottlenecked populations.
# Here, I am extracting all statistics.
# In this script: low genetic variation, initially large, density independent.
# SN - adapted from simulations/run_final_simulations/sim_alldata_n100_a000_lovar.R
# adapted and run on September 22 2020, re-run January 2021,
# archived September 29, 2022

### Clear namespace
rm(list = ls())

### Load packages

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

### Load source materials

# Get simulation functions
source('sim_functions.R')

# Define number of trials
trials = 4000

# Define parameters
pars = data.frame(
  trial = 1:trials,
  n.pop0 = 100,
  end.time = 15,
  init.row = 1e4,
  n.loci = 25,
  w.max = 2,
  theta = 2.75,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0
)

liszt = vector('list', nrow(pars))

set.seed(541002)

for (i in 1:nrow(pars)) {
  sim.output = sim( a = c(1/2, -1/2), params = pars[i,] )
  
  demo.summ = sim.output %>%
    group_by(gen) %>%
    summarise(n = n(),
              gbar = mean(g_i),
              zbar = mean(z_i),
              wbar = mean(w_i))
  
  gene.summ = sim.output %>%
    select(-c(g_i, z_i, w_i, r_i, fem)) %>%
    gather(key = loc.copy, value = val, -c(i, gen)) %>%
    mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
    group_by(gen, locus) %>%
    summarise(p = mean(val > 0)) %>%
    group_by(gen) %>%
    summarise(p.fix.pos = mean(p == 1),
              p.fix.neg = mean(p == 0),
              v = sum(2 * p * (1 - p)) / pars$n.loci[1])
  
  liszt[[i]] = cbind(demo.summ, gene.summ %>% select(-gen)) 
  
  print(paste0('did 100 hi var ', i, ' of ', nrow(pars)))
}

all.sims = unroller(liszt) %>%
  merge(y = pars %>% select(trial, n.pop0, alpha) %>% mutate(low.var = FALSE), by = 'trial') %>%
  group_by(trial)  %>%
  mutate(ext.gen = max(gen),
         extinct = !ext.gen %in% pars$end.time[1]) %>%
  ungroup()

write.csv(all.sims, row.names = FALSE,
          file = "simulations/outputs/alldata_n100_a000_hivar.csv")

sessionInfo()

# session info 29 Sept. 2022

# R version 4.1.2 (2021-11-01)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.6 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] tidyselect_1.1.1 tidyr_1.1.3      dplyr_1.0.7      ggplot2_3.3.5   
# 
# loaded via a namespace (and not attached):
#   [1] fansi_0.5.0      withr_2.4.2      crayon_1.4.1     utf8_1.2.2       grid_4.1.2       R6_2.5.0         lifecycle_1.0.0  gtable_0.3.0     magrittr_2.0.1  
# [10] scales_1.1.1     pillar_1.6.2     rlang_0.4.11     generics_0.1.0   vctrs_0.3.8      ellipsis_0.3.2   tools_4.1.2      glue_1.4.2       purrr_0.3.4     
# [19] munsell_0.5.0    compiler_4.1.2   pkgconfig_2.0.3  colorspace_2.0-2 tibble_3.1.3  