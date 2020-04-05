library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

n.loci = 40

load('~/demgen/model_building_eg_n0_20.RData')

unroller = function(sim.list) {
    sim.list %>%
        do.call(what = rbind) %>%
        mutate(trial = cumsum(c(1, as.numeric(diff(gen) < 0))))
}

all.sims = unroller(liszt) %>%
    select(trial, names(.)[!grepl('trial', names(.))])

all.sims %>%
    select(-c(w_i, r_i, fem)) %>%
    gather(key = allele.copy, value = gen.val, -c(trial, i, b_i, gen)) %>%
    mutate(allele = allele.copy %>% gsub(pattern = '^[ab]', replacement = '') %>% as.numeric()) %>%
    group_by(trial, gen, allele) %>%
    summarise(p.pos = mean(gen.val > 0)) %>%
    mutate(sigma.a = 2 * p.pos * (1 - p.pos) / n.loci) %>%
    group_by(trial, gen) %>%
    summarise(sigma.a = sum(sigma.a)) %>%
    write.csv('~/demgen/allelic_all_sigma_a.csv', row.names = FALSE)


