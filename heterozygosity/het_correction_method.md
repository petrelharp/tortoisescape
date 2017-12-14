# Heterozygosity

We would like to estimate the per-individual heterozygosity, given

1. Proportion of mismatching, overlapping sites on different reads, $x$.
2. Proportion of mismatching, overlapping sites on the same read, $z$.
3. Quality score distributions.

If the per-base error rate is $\epsilon$ and the true heterozygosity is $h$,
then $x = 2 \epsilon (1-\epsilon) + ((1-\epsilon)^2 + \epsilon^2) h/2$,
so given an estimate of $\epsilon$, we will estimate
$$ h = 2 \frac{(x - 2 \epsilon (1-\epsilon))}{(1-\epsilon)^2 + \epsilon^2}. $$

## Error rate estimates

Suppose that quality scores are *relatively* well calibrated,
i.e., that a base with quality score $Q$ in a given individual has probability
$$ a \times 10^{-Q/10}, $$
of being in error, for some constant $a$.
Also suppose that among overlapping sites, the proportion of sites with quality $Q$ is $p(Q)$.
Then the overall density of errors (ignoring the contribution of double errors) is
$$
    z = \sum_Q p(Q) a \times 10^{-Q/10} .
$$
Based on this, we would estimate
$$
    a = \frac{z}{\sum_Q p(Q) 10^{-Q/10}} .
$$

Based on this, we want to now predict the error rates for the sites we actually used in computing $x$ above.
Suppose that the mean read length, after merging overlapping reads, is $L$,
and the mean overlap length is $K$, so a proportion $K/L$ of our bases are produced with two overlapping reads, 
while $1-K/L$ are based only a single read.
Suppose that among the non-overlapping reads, the proportion of sites with quality $Q$ is $q(Q)$.
Then we predict the mean per-base error rate to be
$$
    \epsilon = (1-K/L) \sum_Q q(Q) a \times 10^{-Q/10} ,
$$
again ignoring any error in the bases produced using two overlapping reads.


