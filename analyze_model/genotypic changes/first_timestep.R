# Generate 2000 simulation runs on server.
# Do only first tim step; I want to see the allelic change.
# Run May 11, 2020

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

rm(list = ls())

setwd('~/demall')

source('base_model_source/sim_functions.R')

n.trials = 2000

liszt = vector(mode = 'list', length = n.trials)

set.seed(603185)

for (trial in 1:n.trials) {
  liszt[[trial]] = sim(
    a = c(-1/2, 1/2),
    params = data.frame(end.time = 2,
                        init.row = 1e4,
                        n.loci = 25, 
                        n.pop0 = 20,
                        w.max = 2, 
                        theta = 2.6, 
                        wfitn = sqrt(1 / 0.14 / 2),
                        sig.e = 0)
  ) %>%
    mutate(trial = trial)
}


liszt %>%
  do.call(what = rbind) %>%
  write.csv('out/first_timers.csv',
            row.names = FALSE)

