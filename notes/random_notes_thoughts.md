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


### Density dependence as a variable

(Nov. 18)

It's worth mentioning in the discussion that the strength of density dependence (alpha) is itself a variable. Increasing the resource or habitat availability will reduce the strength of density dependence and avert some of its negative consequences for rescue. Also as Mike K. suggested, perhaps the strength of density dependence can be thought of as a heritable trait, itself subject to evolution.

### Overlapping generations

(Dec. 28)

Been thinking... so many models of rescue and adaptation assume there are no overlapping generations. This makes modeling phenotypic change much simpler. But we know that lots of organisms are perennials with overlapping generations. It would be very easy to flexibly adapt our model to compare discrete generations with overlapping generations.

The easiest way to do this in simulations is to make the beetles biennial; here, cut the *annual* fitness of females in half but have them survive for an extra year. Their per-capita lifetime fitness should be the same... is it necessarily true that the population growth rate is the same though? Think harder about this. (Individuals age 1 contribute 1 + E[W_i] to next generation because they survive... - does it make sense to have survival proportional to fitness?). It would be a little difficult to  have a longer lifespan with stochastic survival, as there is some finagling to set growth rates to be equivalent to the annual life history model.

One note on modeling biennials - a lot of plants are biennial, including some weedy species. This could be interesting for pest control.

I haven't looked too seriously into what is known about the rate of phenotypic change with overlapping generations. A. Prugel-Bennett has some models (and analytical work) demonstrating that adaptation is slowed down. This is pretty intuitive - maladapted alleles stick around for longer and will mix with adapted alleles, meaning maladapted alleles remain in the population for longer. I should look into some textbooks (Lynch and Walsh, Falconer and Mackay) to see what they say. Can the rate of adaptation be analytically characterized in biennials? 

A single sudden environmental shift would be trivial (likely) - adaptation would be slower and extinction is probably more likely. To make this project more interesting, it would be interesting to focus on a stochastic environment. Stochasticity is potentially what is making longer survival beneficial. Then, we may end up in cases similar to the Lyberger et al. paper where higher heritability (or greater levels of phenotypic change) increase extinction risk by having the population adapt to noise. This probably makes sense to think about a moving target with stochasticity, i.e., $\theta_t = \theta_0 + (b)(t) + e_t$ for random variable $e_t$. Here there's likely some classic lit to refer to about survival strategies in varying environments. Nancy and Dan would be great to talk to here.

##### Intrinsic fitness and population growth

Thing that immediately comes to mind is that with overlapping generations, intrinsic fitness is decoupled (somewhat) from population growth. Say an individual lives for exactly two years (i.e., survival probability after year 1 is unity) and has fitness W_i. Its expected lifetime fitness is (2*W_i). This is the same as in an annual model where R_i ~ Pois(2*W_i). 

However, the population growth rate is different. Fitness is how many copies of the genotype are produced over the course of one's lifetime. Population growth rate (per capita) is how many individuals each individual contributes to the next generation. If survival probability is fixed at 1, then individuals contribute $1 + E[W_i]$ to each generation. Thus, if $p_1$ is the proportion of individuals at age 1 and there are no fitness differences between generations, then the growth rate should be something like $p_1 (1 + W_1) + (1 - p_1) W_2$. It seems like a reasonable assumption that $W_2 > W_1$. If population growth is negative than it also seems logical that $p_1 > p_2$ and the opposite if the population is growing. 

(Also is it true that, for X ~ Pois(L), Y ~ Pois(2L), that E[X1 + X2] = E[Y]? It may not be...)

It's worth wondering what happens if fitness also determines survival probabilities. If $K_i = exp(-(z_i)^2/(2*w^2)) = W_i / W_{max}$, then the expected contribution of a 1 year old to the next generation is $K_i + W_i = W_{max} (K_i + 1)$. The expected lifetime fitness is $W_i (1 + K_i)$. It would be interesting to model these as two separate traits.


