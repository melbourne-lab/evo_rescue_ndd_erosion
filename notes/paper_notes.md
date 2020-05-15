### Reznick, 2016. Hard and soft selection revisited: How evolution by natural selection works in the real world. Journal of Heredity, 2016, 3 - 14.

##### Haldane and the "classical hypothesis"
- Models assume there is one best genotype (w/ highest fitness)
	- difference between this fitness and (mean) observed fitness is *genetic load*
- Assume that adapting to optimum occurs through slow allelic replacement (one at a time)
	- this is the "classical hypothesis"
	- somehow (?) it follows from this that genetic variation is limited (presumably because it is selected against)
		- although, we know that there is large natural variation

##### Hard versus soft selecton
- Haldane's selection is hard
- Soft selection: environment has vacancies to be filled, w/ genotypes competing to fill them; it will be filled either way
	- this somehow could explain genetic variation...?
	- in a sense these models assume population growth is exponential and ignore resource limitations

##### Guppy example
- Naturally occurring high predation (HP) and low predation (LP) environments
- Predation risk higher for larger individuals
- In HP environments, smaller reproductive individuals, more reproduction
- However, no evidence of age-dependent mortality changing across predator-environments (i.e., no age x predator interaction on mortality)
	- the population ecology (size distribution) explains this: in LP environments, resource stress induced by unproliferated growth not present in HP environments
- Bassar et al., 2010, 2012, 2013, experimentally found that HP guppies have higher fitness than LP at low population density, but not at high density
	- i.e., adaptation to density influences adapting to LP environments, because guppies in LP environments are less sensitive to density

### Shaw and Shaw, 2014. Quantitative genetic study of the adaptive process. Heredity, 112: 13 - 20.

##### Fisher's Fundamental Theorem of Natural Selection
- Change in mean fitness is additive genetic variance divided by mean absolute fitness
	- V_A(W) / \bar{W}
	- doesn't account for mutations (short time scale)
- Idea hasn't been used empirically very much
- Charlesworth, 1987, suggests VA should be slight at equilibrium

##### Quantitative genetic models
- Turelli and Barton, 1984 demonstrate that with strong selection, dist'n of breeding values is approximately normal (we should confirm)
- Under Gaussian selection (like our models), phenotype additive environmental and genetic (e + a = z)
	- Distribution of genotypic effects on absolute fitness, E(W|a)
	- check appendix for PDFs of W?
- "the inclusion in the model of environmental variance as a contributor to variation in the trait causes the genotypic distribution for fitness to be shifted toward intermediate values. In that case, no genotype has an expected fitness that coincides with the maximal fitness, and genotypes with low fitness values in the absence of higher variation have higher expected fitness when it is included"
	- so environmental variance speeds up when far away but slows down when you get closer (because it's adding noise)
- I'm really lost by the second paragraph of the "Fitness under Gaussian selection" paragraph...
	- What is V_A(W) and what is V_G(W)?
		- V_G(W) is the variance of the PDF of W, f(W)
		- V_A(W) is the additive genetic variance for fitness ?
	- Fig 2: increasing sig_a (additive genetic variance) increases V_G(W)
- Pease et al., 1989, has a population experiencing gaussian selection with a moving target
- Burger and Lynch, 1995, look at persistence of populations using environmental change
- Ronce et al., 2009, has a model very similar to ours, as well as a framework for incorporating inbreeding
	- this model allows genetic variance to evolve freely
	- includes mutations
	- here, a rapid increase in trait V_G(W) with a subsequent change in sig^2_a... is this due to the moving target? What is going on here.
	- also calculate V_A(W) somehow... but how?
	- take home of this is that V_A(W) can increase when the selection regime changes

### Ronce, Shaw, Rousset, and Shaw. Is inbreeding depression lower in maladapted populations? A quantitative genetics model. Evolution, 63(7): 2009.

##### How do ecological scenarios and influence inbreeding depression?
- A widely-accepted model of inbreeding: two alleles, one is unconditionally deleterious, at least partially recessive
	- see Charlesworth and Charlesworth, 1999
	- Acc'd to this model, inbreding should purge deleterious aleles
- These models directly link genotype to fitness, no traits, unlike qGen
- This paper looks at inbreeding w/ stabilizing selection on quantitative traits
	- see Lande and Schemske, 1985, for something simpler

##### Model
- trait, z, has additive action of n loci, diploid, panmictic
	- z = mu + sum_i=1^n (am_i + af_i) + e
	- am_i, af_i are allelic effects in herited from parents, mu is mean
- trait under Gaussian selection
- I = 1 - (\bar{W}_s) / (\bar{W}_r)
	- I = inbreeding depression
	- selfing versus random mating fitesses (selfing: self is parents?)
- Analytic model follows Lande, 1976
	- here, an expression for inbreeding depression! interesting

Come back to this paper later. This is going in some directions beyond what we need right now.

### Osmond and de Mazancourt, 2013. How competition affects evolutionary rescue. Phil Trans R Soc B 368: 20120085.

##### Competition
- Reduces population abundance (extinction more likely)
- But, also increase selective pressure, increase rate of evolution (refs herein)
- What is the influence of competition (intra and interspecific) on rescue?

##### One species model
- Asexual population: density- and frequency-dependent growth
- Continuous traits? Continuous time
- Growth of # individuals w/ a trait is logistic
	- features per capita intrinsic growth rate
	- carrying capacity for each trait is a defined by gaussian
- Intraspecific competition: same trait value compete most strongly
- Mutation is rare, s/t demography is at equilibrium, evolution is slow
- Population initialized at constant trait value, perfectly adapted, reaching optimal carrying capacity
	- environment changes, carrying capacity changes
- N_c is extinction threshold
- for a given N_c, get an expression for expected time (t_R) to surpass N_c
	- t_R increases with increasing environmental change and N_c
	- unimodial w/r/t environmental tolerance (sigma_K)/selection strength
		- tol too high, and population is never at extinction risk
		- tol too low and population immediately extinct
- One term in here is the analogue to the rate of change in the trait value from G&H
	- in O&dM model, evolution stops with absence of mutations

##### Discussion
- Competition which increases the strength of selection may speed up adaptation... where is this coming from? Ref Schluter et al., 1985.
- Refs: "intraspecfic competition often has relatively little impact on selectiv pressures" (refs)
- Impacts of NDD:
	- Boulding and Hay: populations stay small, where NDD doesn't act.
	- Chevin and Lande: NDD *makes* populations small, where extinction is more likely
- This model assumes "evolution is slow, and hence, on the time-scale of evolution, populations are always at carrying capacity"
	- this paper is different from our study b/c we're interested in *rapid adaptation* and the stochasticity is important
	- "we assume greater abundances lead to faster evolution because they lead to greater mutational inputs" - population size is important for mutation
	- "intrinsic growth rate is therefore irrelevant for evolutionary rescue in our mode" (due to the demographic equilibrium)
		- this is evolutionary rescue under K-selection rather than r-selection

So... this paper is all about how mutations drive rescue and how that interacts with density dependence. Model is fundamentally different because populations are assumed to be at carrying capacity for most of the time... but then how to explain Fig. 2? Not sure.

### Vinton, A.C., and Vasseur, D.A. Evolutionary rescue is determined by differntial selection on demographic rates and density dependence. bioRxiv preprint.

##### Evolutionary rescue
- Four factors influencing potential for rescue (Bell and Gonzalez, 2009)
	- Initial population size (N_0)
	- Genetic variability (standing variation and mutations)
	- genetic variability due to dispersal (interesting...)
	- extent/severity of environmental change
	- (n.b., three of these are in G&H's equations)

##### Demography
- Models typically only focus on $r$, ignore separate births nad deaths
- e.g., Holt, 1990: populations with high birth and death vs. low birth and death
	- these have the same $r$, but different rates of turnover, mutation, etc.

##### Model
- Continuous time, ndividual-based logistic growth
- N' / n = b(gamma0, gamma1) - d(gamma0, gamma1) where
	- gamma0: density independent contributions
	- gamma1: density dependent contributions
	- b(), d() depend on rates, carrying capacity K
	- carrying capacity assumes all individuals well adapted
- environmental effect: mu_0(t,f), varies sinusoidally over t
	- e_mu = abs(mu - mu_0) is a measure of maladaptation
- Birth and death models:
	- Case 1a: environment alters b(), density independent modification
	- Case 1b: environment alters b(), density dependent modification
	- i.e., the degree to which b() is lowered does/does not depend on environment
	- Case 2a, 2b: density modifies d(), density independent/dependent
	- Density dependence: only see decreases/increases in rates when pop'n size is large
		- e.g., little maladaptation at size one
- Stochastic implementation via the Gillespie algorithm as birth/death process
	- good description of Gillespie here!
	- includes mutation, randomly drawn from U(-0.3,0.3), added to parential trait

##### Results
- r is differece of birth and death, but stochasticity is determined by sum (?)
	- i.e., variance is due to the sum of variances in birth and death?
	- refs herein
- death models have higher demographic stochasticity
	- there's a reason here but I don't understand it.
- due to model construction, density dependence means *less* variable population size
	- (by construction this means you maintain maladapted individuals at low size)
	- notably, this is is relative to a DI-*effect* on rates, not on DI itself
	- this is a case with phenotype x growth rate interaction
- increasing persistence time with a faster-changing environment?

(didn't quite finish, but worth re-reading)

