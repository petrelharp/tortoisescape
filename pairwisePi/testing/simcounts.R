#!/usr/bin/Rscript

pwp <- function (counts,coverages) {
    pwp <- numeric( ncol(coverages)^2 )
    dim(pwp) <- c( ncol(coverages), ncol(coverages) )
    for (i in 1:ncol(coverages)) {
        for (j in 1:ncol(coverages)) {
            weights <- sqrt(coverages[,i]) * sqrt(coverages[,j])
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
# binomial errors

nloci <- 1e6

mean.coverage <- 2
phet <- .008  # per site
perr <- .001  # per read

hets <- cbind( rbinom( nloci, size=1, prob=phet ), rbinom( nloci, size=1, prob=phet ) )
coverages <- do.call( cbind, lapply( 1:2, function (k) {
            x <- rpois( nloci*2, mean.coverage )
            x[x>0][1:nloci]
        } ) )
counts <- sapply( 1:2, function (k) {
        rbinom( nloci, size=coverages[,k], prob=ifelse(hets[,k],1/2,perr) )
    } )

all.pwp <- pwp(counts,coverages)
coversites <- lapply( 1:6, function (k) {
        coverages[,1]==coverages[,2] & coverages[,1]==k
    } )
pwp.list <- lapply( coversites, function (kk) {
        pwp(counts[kk,],coverages[kk,])
    } )
sapply( coversites, sum )


####
# Beta-binomial errors

nloci <- 1e6

mean.coverage <- 2
phet <- .008  # per site
vhet <- 4*phet # controls correlation of heterozygous sites
perr <- .001  # per read
verr <- 10

hets <- rbinom( nloci, size=1, prob=phet )
hets <- cbind( hets, rbinom( nloci, size=1, prob=(phet+hets)/(phet+1) ) )
# perrs <- cbind( rbeta( nloci, perr*verr, (1-perr)*verr ), rbeta( nloci, perr*verr, (1-perr)*verr ) )  # SEPARATE error probs per indiv
perrs <- rbeta( nloci, perr*verr, (1-perr)*verr ); perrs <- cbind( perrs, perrs )  # SAME error probs per indiv
coverages <- do.call( cbind, lapply( 1:2, function (k) {
            x <- rpois( nloci*2, mean.coverage )
            x[x>0][1:nloci]
        } ) )
counts <- sapply( 1:2, function (k) {
        rbinom( nloci, size=coverages[,k], prob=ifelse(hets[,k],1/2,perrs[,k]) )
    } )

all.pwp <- pwp(counts,coverages)
coversites <- lapply( 1:4, function (k) {
        coverages[,1]==coverages[,2] & coverages[,1]==k
    } )
pwp.list <- lapply( coversites, function (kk) {
        pwp(counts[kk,],coverages[kk,])
    } )
sapply( coversites, sum )

