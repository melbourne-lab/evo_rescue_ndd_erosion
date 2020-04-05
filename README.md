# Allele-specific models of demographic and genetic rescue

This model is an allele-specific model of evolutionary rescue. Previously, we had been using the quantitative-trait model of Gomulkiewicz and Holt, 1995, for modeling a population adapting to novel conditions. However, this model was not exactly compatible with the genetics of our system, as the weighting of breeding values due to selection on fecundity produced a genetic variance that was not constant over time.

This led to the question, how should genetic variance change over time? One would expect that in the absence of mutation, selection would reduce variance as less fit genotypes disappear from the population. The rate of this disappearance, and the disappearance itself, are hard to characterize.

To investigate this I have constructed an allele-specific model of evolving populations of flour beetles. The model features m loci, which combine additively to produce a genotype, which is then used to produce a phenotype. Contributions of each locus are scaled such that, in the limit of large m, breeding values are initially normally distributd (according to the Central Limit Theorem). When individuals mate, alleles are randomly segregated to offspring. Full model write-up coming soon.

### Folder contents

Current repo structure

- `init_model` directory containing files from initial versions of model constructed in March, 2020
