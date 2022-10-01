# 15-generation simulations of rescue, recording all population-level variables
# In this script: high genetic variation, initially small, density independent.
# adapted from simulations/run_final_simulations/sim_alldata_n100_a000_lovar.R
# first adapted and run on September 22 2020, re-run January 25 2021,
# re-run for archiving October 1 2022

### Clear namespace
rm(list = ls())

### Load packages
library(ggplot2)
library(dplyr)
library(tidyr)

### Load source materials

# Get simulation functions
source('sim_functions.R')

# Define number of trials
trials = 4000

# Define parameters
pars = data.frame(
  trial = 1:trials,
  n.pop0 = 20,
  end.time = 16,
  init.row = 1e4,
  n.loci = 25,
  w.max = 2,
  theta = 2.8,
  wfitn = sqrt(3.5),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0
)

# Initialize list for storage
liszt = vector('list', nrow(pars))

### Run simulations
set.seed(707106)

for (i in 1:nrow(pars)) {
  
  sim.output = sim( params = pars[i,] )
  
  # Get demographic variables
  demo.summ = sim.output %>%
    group_by(gen) %>%
    summarise(
      n = n(),
      gbar = mean(g_i),
      zbar = mean(z_i),
      wbar = mean(w_i),
      pfem = mean(fem)
    )
  
  # Get genetic variables
  gene.summ = sim.output %>%
    select(-c(g_i, z_i, w_i, r_i, fem)) %>%
    gather(key = loc.copy, value = val, -c(i, gen)) %>%
    mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
    group_by(gen, locus) %>%
    summarise(p = mean(val > 0)) %>%
    group_by(gen) %>%
    summarise(
      p.fix.pos = mean(p == 1),
      p.fix.neg = mean(p == 0),
      v = sum(2 * p * (1 - p)) / pars$n.loci[1]
    )
  
  liszt[[i]] = cbind(demo.summ, gene.summ %>% select(-gen)) %>%
    mutate(trial = i)
  
  print(paste0('did 20 hi var ', i, ' of ', nrow(pars)))
}

### Process and export

# Convert list to data frame and add parameters
all.sims = liszt %>%
  do.call(what = rbind) %>%
  merge(y = pars %>% select(trial, n.pop0, alpha) %>% mutate(low.var = FALSE), by = 'trial') %>%
  group_by(trial)  %>%
  mutate(
    ext.gen = max(gen),
    extinct = !ext.gen %in% pars$end.time[1]
  ) %>%
  ungroup()

# Export CSV
write.csv(all.sims, row.names = FALSE,
          file = "simulations/outputs/alldata_n20_a000_hivar.csv")

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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] tidyr_1.1.3   dplyr_1.0.7   ggplot2_3.3.5
# 
# loaded via a namespace (and not attached):
#   [1] fansi_0.5.0      withr_2.4.2      crayon_1.4.1     utf8_1.2.2      
# [5] grid_4.1.2       R6_2.5.0         lifecycle_1.0.0  gtable_0.3.0    
# [9] magrittr_2.0.1   scales_1.1.1     pillar_1.6.2     rlang_0.4.11    
# [13] rstudioapi_0.13  generics_0.1.0   vctrs_0.3.8      ellipsis_0.3.2  
# [17] glue_1.4.2       purrr_0.3.4      munsell_0.5.0    compiler_4.1.2  
# [21] pkgconfig_2.0.3  colorspace_2.0-2 tidyselect_1.1.1 tibble_3.1.3    