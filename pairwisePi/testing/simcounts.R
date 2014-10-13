#!/usr/bin/Rscript

pwp <- function (counts,coverages) {
    pwp <- numeric( ncol(coverages)^2 )
    dim(pwp) <- c( ncol(coverages), ncol(coverages) )
    for (i in 1:ncol(coverages)) {
        for (j in 1:ncol(coverages)) {
            weights <- (coverages[,i]) * (coverages[,j] - if (i==j){1}else{0} )
            if (i!=j) {
                probs <- counts[,i]/coverages[,i] * (1-counts[,j]/coverages[,j]) + (1-counts[,i]/coverages[,i]) * counts[,j]/coverages[,j]
                pwp[i,j] <- sum( weights*probs, na.rm=TRUE ) / sum(weights)
            } else {
                usethese <- ( coverages[,i] > 1 )
                probs <- 2 * counts[,i] * (coverages[,i]-counts[,i]) / ( coverages[,i] * (coverages[,i]-1) )
                pwp[i,i] <- sum( (weights*probs)[usethese] ) / sum(weights[usethese])
            }
        }
    }
    return(pwp)
}

####
# simple; binomial errors, same set of heterozygous sites in all samples

nloci <- 1e6
nsamples <- 4

mean.coverage <- 2
phet <- .08  # per site
perr <- .01  # per read

# same segregating sites between torts
hets <- rbinom( nloci, size=1, prob=phet )
coverages <- do.call( cbind, lapply( 1:nsamples, function (k) { 
            rpois( nloci, mean.coverage ) 
        } ) )
counts <- sapply( 1:nsamples, function (k) {
        rbinom( nloci, size=coverages[,k], prob=ifelse(hets,1/2,perr) )
    } )

# right answer should be
true.pwp <- matrix( phet * 1/2 + (1-phet) * 2 * perr * (1-perr) )
# and we get
all.pwp <- pwp(counts,coverages)
# yep!
all.pwp / as.numeric(true.pwp)

coversites <- lapply( 1:6, function (k) {
        coverages[,1]==coverages[,2] & coverages[,1]==k
    } )
pwp.list <- lapply( coversites, function (kk) {
        pwp(counts[kk,],coverages[kk,])
    } )

# yep!
rbind( sapply( pwp.list, function (x) range( x/all.pwp, na.rm=TRUE ) ),
    sample.size=sapply( coversites, sum ) )



####
# Correlated heterozygous sites and correlated error probabilities

nloci <- 1e6
nsamples <- 4

mean.coverage <- 2
phet <- .08  # per site
vhet <- 4*phet # controls correlation of heterozygous sites
perr <- .01  # per read
verr <- 10   # correlation of errors

hetprobs <- rbeta( nloci, phet*vhet, (1-phet)*vhet )
hets <- do.call( cbind, lapply( 1:nsamples, function (k) {
            rbinom( nloci, size=1, prob=hetprobs )
        } ) )
errprobs <- rbeta( nloci, perr*verr, (1-perr)*verr )

coverages <- do.call( cbind, lapply( 1:nsamples, function (k) {
            rpois( nloci, mean.coverage )
        } ) )

counts <- sapply( 1:nsamples, function (k) {
        rbinom( nloci, size=coverages[,k], prob=ifelse(hets[,k],1/2,errprobs) )
    } )

# our guess
all.pwp <- pwp(counts,coverages)
# correct answers
true.self.pwp <- phet * 1/2 + (1-phet) * 2 * perr * ((1-perr)*verr)/(verr+1)
true.betw.pwp <- (
        phet * (phet*vhet+1)/(vhet+1) * 1/2  +      # both are het
        2 * phet * ((1-phet)*vhet)/(vhet+1) * 1/2 + # one is het, other is not
        (1-phet) * ((1-phet)*vhet+1)/(vhet+1) * 2 * perr * ((1-perr)*verr)/(verr+1) # neither
    )
true.pwp <- matrix( true.betw.pwp, nrow=nsamples, ncol=nsamples ); diag(true.pwp) <- true.self.pwp

# yep!
all.pwp/true.pwp

coversites <- lapply( 1:4, function (k) {
        coverages[,1]==coverages[,2] & coverages[,1]==k
    } )
pwp.list <- lapply( coversites, function (kk) {
        pwp(counts[kk,],coverages[kk,])
    } )

# yep!
rbind( sapply( pwp.list, function (x) range( x/all.pwp, na.rm=TRUE ) ),
    sample.size=sapply( coversites, sum ) )


####
# Correlated heterozygous sites and correlated error probabilities
#  also with mapping ambiguity

nloci <- 1e6
nsamples <- 4

mean.coverage <- sort( 6*runif(nsamples) )
phet <- .08  # per site
vhet <- 4*phet # controls correlation of heterozygous sites
perr <- .01  # per read
verr <- 10   # correlation of errors
pmap <- .2   # contract this proportion of sites

hetprobs <- rbeta( nloci, phet*vhet, (1-phet)*vhet )
hets <- do.call( cbind, lapply( 1:nsamples, function (k) {
            rbinom( nloci, size=1, prob=hetprobs )
        } ) )
errprobs <- rbeta( nloci, perr*verr, (1-perr)*verr )

orig.coverages <- do.call( cbind, lapply( 1:nsamples, function (k) {
            rpois( nloci, mean.coverage[k] )
        } ) )

# randomize over alleles
polarization <- rbinom( nloci, size=1, prob=1/2 )
orig.counts <- sapply( 1:nsamples, function (k) {
        pp <- ifelse(hets[,k],1/2,errprobs)
        rbinom( nloci, size=orig.coverages[,k], prob=ifelse(polarization,pp,1-pp) )
    } )

# sum these sites with preceding good site
badmap <- c(0, rbinom( nloci-1, size=1, prob=pmap ) )
counts <- sapply( 1:nsamples, function (k) {
        diff( c(0,cumsum( orig.counts[,k] )[badmap==0] ) )
    } )
coverages <- sapply( 1:nsamples, function (k) {
        diff( c(0,cumsum( orig.coverages[,k] )[badmap==0] ) )
    } )

# see, het and coverage correlated
tapply(counts[,1]>0,coverages[,1],mean)

# our guess
all.pwp <- pwp(counts,coverages)

# yep! it does NOT depend on coverage of the sample
all.pwp

## NOT YET:
# correct answers
true.map <- sum( (1-pmap) * pmap^(0:20) * 1/(1:21) )  # prob unifly chosen pair of reads at a the same site (not uniform across sites!) map to same spot

true.self.pwp <- ( 
        true.map * ( phet * 1/2 + (1-phet) * 2 * perr * ((1-perr)*verr)/(verr+1) ) + # both reads map to same spot
        (1-true.map) * 1/2  # otherwise, completely independent
    )

# copied from above
true.betw.pwp <- (
        phet * (phet*vhet+1)/(vhet+1) * 1/2  +      # both are het
        2 * phet * ((1-phet)*vhet)/(vhet+1) * 1/2 + # one is het, other is not
        (1-phet) * ((1-phet)*vhet+1)/(vhet+1) * 2 * perr * ((1-perr)*verr)/(verr+1) # neither
    )
true.pwp <- matrix( true.betw.pwp, nrow=nsamples, ncol=nsamples ); diag(true.pwp) <- true.self.pwp

# not yet on the math
all.pwp/true.pwp

