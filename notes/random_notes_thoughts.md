# Random demgen notes

Keeping stray notes relevant to the allele-specific demographic rescue project in this document.

#### Within-locus and among-locus diversity/variance

(June 9)

Was thinking today about why within-locus variance in allele frequencies (i.e., allele frequencies closest to 0.5) is important for changes in population mean genotype. This relates to Lande's (1976) derivation showing that allelic change increases with phenotypic variance. This problem makes more sense when you think about why allelic change may be *lowest* while allele frequencies are close to 0 or 1. In this case, there are very few copies of one allele which can self-replicate (through sexual reproduction) to increase their own frequency. This argument needs a little bit more rigor (as it says nothing about the magnitudes of changes in allele frequencies *towards* 0 or 1), but is a good illustration.

This brings up a good point - there are two separate genetic processes here governing the changes in allele frequencies! One is a stochastic drift-like process outlined above which depends on the allele frequencies and the other is a directional shift or stabilization induced by selection. The relative sizes of these do not necessarily correlate with each other. E.g., you can have a case where the positive allele is favored but at very low frequency; here, the effect of selection is strong (the selection strength is strong because the population mean is further from the phenotypic optimum) but the within-locus variance is small (one allele dominates) slowing the progress of allele frequency change.

Why is the loss of variance at one locus (i.e., fixation) problematic? Overall phenotypic change is the sum of changes occurring at individual loci. Because of directional selection, the mean (expected) phenotypic change here is positive, even if selection means there may be some negative changes (i.e., loci with allele frequencies moving away from the optimum). Fixation (variance of zero) means that a locus is incapable of contributing to this mean, dragging the sum of changes closer to zero.

#### Simulation the non-evolving treatment

(June 9)

Was thinking about this earlier today. If we are trying to simulate all four of our experimental treatments, how can we simulate a population which is not evolving? The key here is that we want *only* demographic stochasticity, i.e., population size. My sense is that if we specify the initial allele frequencies (default is p = 0.5 for each locus) that for each generation, we should just randomly reassign the allele frequencies to match with the initial allele frequencies.

#### Environmentally-induced phenotypic variance

(June 24)

A question. Thinking about how to incorporate bottlenecks into simulated populations. Perhaps this question has a trivial answer, or thinking about it here will inspire me. If we incorporate a poulation from an old environment into a new environment, does the phenotype from the old environment necessarily cary over? Or will there be a new trait value, determined in part with noise from the new environment? 


