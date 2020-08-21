## Introduction

#### Evolutionary rescue is...
	- Population adapted to a set of conditions
	- Environment changes, now mean fitness below replacement rate
	- Population will go extinct (deterministically) unless population adapts to the point that lambda > 1
	- This is relevant considering widespread anthropogenic climate change and habitat destruction/degradation

#### Theory of evolutionary rescue: we know what is important through a theory making a lot of simplifying assumptions
	- Gomulkiewicz and Holt (1995) model
		- Classic "u" shaped trajectory
		- Determinants of extinction risk: initial population size, amount of genetic variance, heritability, degree of maladaptation
		- Model assumptions: quantitative trait, density independence
	- Unckless, Orr, others
		- Assumptions about genomic architecture (e.g., ploidy, complication of trait)
	- Other models (e.g., Barfield and Holt, 2016)

#### Evolutionary rescue is of particular importance for conservation targets
	- Features of cons. targets include: ploidy, sexual reproduction, more complicated demography
	- Additionally, models are often deterministic, or lacking many realistic elements of stochasticity such as demographic stochasticity from small population size and genetic stochasticity through the random segregation of alleles and combinations of multiple traits
	- We can make better predictions and forecasts and have better management strategies if we are operating under models that are more realistic to the biology of conservation targets

#### Density dependence is important for population dynamics
	- Negative density dependence: population growth is depressed as populations grow large
		- This is typically induced by interspecific competition for resources (e.g., ...), including physical space, alhough it may result from other processes such as cannibalism (cannibalism can likely be left out)
	- Depending on the strength of density dependence, NDD keeps populations at a small-to-medium size for longer
	- This typically is left out of models of rescue and adaptation for tractability and also because by definitions smaller populations experience density dependence less than larger ones
	- However, by keeping populations at smaller size, NDD may increase the risk of extinction
	- Density dependence does not necessarily only affect demographic processes; it may also influence evolution
		- By keeping populations smaller for longer, NDD may increase the probability of inbreeding, genetic drift, and the erosion of genetic variation
		- This has the potential to create a feedback loop
		- This may require modified quantitative genetics models

#### We consider the effects of density dependence on evolutionary rescue by combining and modifying pre-existing models
	- Specifically, we combine a modified stochastic Ricker model with a finite-locus quantitative genetic model
	- The demographic model
		- Stochastic implementation includes many sources of stochasticity
			- Sex, demographic stoch., inter-individual variance in fecundity, 
			- These are important: see Melbourne and Hastings, 2008
		- Negative density dependence controlled by alpha parametr
	- Genetic model
		- Diploid, bi-allelic model with finite number of independently segregating loci
		- Stochasticity through random segregation of alleles
	- We examine the effects of density dependence on **various variables (e.g., population size, extinction probability, time to rescue…) need to be specific here**
	- We also consider how density dependence interacts with two important variables from G&H: initial size, genetic variance

###### Questions on introduction

- Mention Triboleum?

## Methods

## Results

#### Mean population size at end of simulations
	-  **Figure 1** population sizes over time 
		- unconditional population size, i.e., includes extinctions
	- As with classic rescue models (G&H 1995), reduced variation induces long-term reduction in population size
	- Density dependence, with our parameterization, has larger effects than low genetic variation
		- This means that populations take a *longer* amount of time to reach their initial size
	- Notably, populations with density dependent growth are much smaller at the end of the simulation
		- Long-term extinction risk not seen in density-independent growth (*this could go into discussion*)
	- Finally, we saw evidence of an interaction between NDD and low genetic variance
		- Low genetic variance has larger effects on population size under NDD (*or phrase it the other way around?*)

#### Extinction: how NDD and interactions influence extinction risk
	- **Figure 2** cumulative and/or instantaneous risk of extinction in each time step

#### Adaptation is slower under NDD due to lost genetic variation
	- **Figure 3** genetic variation over time
		- Note that extinct populations are excluded from this average (survivor effects)
	- In all treatments, genetic variation decreases over time as expected
	- Reflecting genetic drift, smaller populations lose more genetic variation over time, esp. in early time steps (when smallest)
	- All treatments: more genetic variation lost under NDD
		- NDD has larger effect (reduction in variance) in larger populations
	- **Figure 4** mean intrinsic fitness (W) or genotype (z) over time
	- Look! NDD slows down adaptation.
	- **Figure 5** fixation probabilities over time
	- Maladaptive fixtion occurring at one point, adaptive fixation happening at other points


###### Questions on results

- General questions: how quantitative to get with results? Serious analysis (effect sizes, etc.) vs. qualitative descriptions

- Figure 4: fitness or genotype?
	- Advantage to fitness: can show where mean fitness R > 1
		- (although, under NDD, R > 1 DNE lambda > 1)
	- Four panels? Have NDD and DID in the same panel
		- also could have corresponding predicted phenotype over time from Lande/G&H in each panel
		- if we include Lande phenotype/fitness then are we seeing relative influences of genetic erosion and demography?
		- (thing to note about Lande/G&H is that IDK if it includes extinctions... I feel like it doesn't)

- Good order for the adaptation stuff? Likely should be broken into multiple paragraphs

## Discussion

General discussion points to include (order, content, etc. not fully thought out)

#### Realistic density dependence has stronger effects than a bottleneck
	- (*effects here of genetic variation are specific to the model we use, which takes a simplistic view of genomic changes during adaptation*)
	- (*nonetheless... our genetic model is cool and good because it gives some mechanistic of the loss of genetic variation over time, producing similar dynamics to Barfield & Holt 2016)*
	- This is important because...

#### Density dependence: the extinction is longer-term than originally anticipated
	- Under the density independent framework, extinction risk is very low at the end of trials
		- Without density dependence constraining growth, populations grow to be very large with effectively zero risk in a stable environment
	- With density dependent growth, populations are held at smaller size, with larger long-term extinction risk
		- With NDD, population growth is *highest* when small, so populations can still rebound from small size
		- However, being held at a smaller size for longer, and having a longer time until reaching original size, increases risk of extinction due to stochasticity
		- *could mention other types of density dependence, e.g., Allee effects, here, although IDK if that's necessary*

#### Modified model reveals feedbacks between population size/dynamics and rate of adaptation
	- Not in previous models, although see (and read, lulz) Barfield and Holt, 2016
	- This is important, as many populations of interest are small, so this feedback loop is important to consider
	- This feedback gives insight and nuance to interpreting effects of N_0 and genetic variance as G&H point out
		- Our model shows that populations that start at a smaller size lose more genetic variation than those that start large, meaning that N_0 small means less adaptive capacity over time
		- However, because NDD acts more on large populations, we see that NDD has larger effects on the genetic erosion in populations where N_0 is large
		- Interactions between NDD and genetic varation: low genetic variation has stronger effects in populations experiencing NDD due to the feedback loop (extinction vortex)

#### (segue from above) Connections to Soule's extinction vortex (Gilpin and Soule?)
	- Describing Gilpin and Soule's vortex
	- Note: this process occurs independtly of NDD, although NDD makes it worse (steeper)
	- In terms of conservation, this tells us that smaller populations potentially facing scarce resources or finite space are at especially large risk of extinction

#### Some mechanistic intuition from the genetic model
	- Loss of genetic variation happens through fixation and changes in Allele frequencies
	- As the number of loci approaches infinity, the changes in allele frequencies have less effect on genetic variance (which is what we see in quantitative trait models)
	- Our models allow fixation both of adaptive and maladaptive alleles
		- Fig. 5 shows some of this...
		- Maladaptive alleles, which place a cap on the adaptive capacity of the population, reach fixation only during a narrow time window when populations are small
			- it's sampling variation that produces this
			- this also demonstrates the importance of initial population size and speed of adaptation - in the G&H model these determine how small the population will be at its lowest point
		- Adaptive alleles reach fixation over a longer period of time
			- fixation of adaptive alleles does not limit the ultimate adaptive capacity (max genotypic value) the population can reach
	- This lost genetic variation may be especially important, with some more complicated and interesting dynamics, if the environment varies and the optimal phenotypic value has large fluctuations or moves in the opposite direction of the initial change (i.e., theta moves back towards zero).

#### Including more realistic features in models... is good for increasing the accuracy and generating predictions for conservation targets and wild or managed populations experiencing evolutionary rescue
