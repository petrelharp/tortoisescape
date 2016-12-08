# Methods


The fact that mitochondrial DNA exists at higher copy number in cells
means that although coverage of the nuclear genome averaged only 1.??x,
coverage of the mitochondria averaged 25??x,
after trimming and filtering reads as described above.
The distribution of coverages across sites in each individual was only very slightly overdispersed
relative to the expected Poisson distribution,
indicating no problems with the reference sequence.
After removing individuals with mean coverage below 15x,
we called genotypes by taking the majority base in each individual at each site.
No sites showed evidence of true heterozygosity (i.e., heteroplasmy),
but some sites had more than one sequenced base in many individuals,
with bases that did not match the individual's genotype showing up in about 4% of the reads
(which is close to the mean nuclear coverage divided by the mean mitochondrial coverage).
This suggests that regions of the nuclear genome are mapping to the mitochondria (i.e., NUMTs, cite),
and indeed,
these sites cluster in long segments along the genome:
Figure [numts] shows the proportion of sequenced bases that do not match the individual's genotype.
Precisely, if $n_{si}$ is the total coverage of site $s$ in individual $i$,
and $k_{si}$ is the number of reads in individual $i$ whose sequenced base at site $s$ 
does not agree with the most commonly sequenced base for that individual,
Figure [numts] shows $\sum_i k_{si} / \sum_i n_{si}$, for each $s$.

![minor-freq-along-mt.pdf](Figure numts: Fraction of sequenced bases along the mitochondrial genome that differ from the genotype.
For each site, the number plotted is $k/n$, where $n$ is the total covarage at that site,
and $k$ is the number of reads whose sequenced base at that site does not agree with the most commonly sequenced base
in that individual.)

Fortunately, NUMTs cover only about 60% of the mitochondrial genome,
so we could use the remaining segments
to estimate the per-sequenced-base error rate remaining after our quality filters
(i.e., the probability that the nucleotide reported at a given site in a given read does not match the true nucleotide).
This is estimated simply as the frequency of sequenced bases not matching the individual's genotype,
$\sum_i k_{si} / \sum_i n_{si}$, averaged across sites not underlying a NUMT.
Across 13,066,260 sequenced nucleotides at 2001 sites in two contiguous stretches, 
only 3016 were putative errors,
giving an estimate of $\epsilon=2.3\times10^{-4}$,
and a standard error of $4.2 \times 10^{-6}$.
There are not enough errors to obtain tight estimates of per-individual error rates,
but per-individual error frequencies are correlated between the two regions
($r=0.55$, $P<10^{-12}$),
and the number of errors per individual are significantly overdispersed relative to the Poisson expectation
(Supplemental figure [error_correlation]).

![error_estimate_correlation.pdf](Figure error_correlation:
The left panel shows the number of erroneous bases against total coverage, per individual.
Red lines give lower and upper 99.9% quantiles based on the Poisson expectation;
the fact that some individuals fall outside of these
provides statistically significant evidence 
that error rates differ significantly between individuals.
The right panel shows error rates in errors per kilobase
estimated separately for each individual in the two distinct regions of the mitochondrial genome
that were combined to get the overall estimate.)


# Results


