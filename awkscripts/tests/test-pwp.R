nsites <- 100
counts <- cbind( x0=rpois(nsites,1), x1=rpois(nsites,1), y0=rpois(nsites,2), y1=rpois(nsites,2), z0=rpois(nsites,10), z1=rpois(nsites,10) )
coverages <- counts[,2*(1:(ncol(counts)/2))] + counts[,2*(1:(ncol(counts)/2))-1]

pwp <- numeric( ncol(coverages)^2 )
dim(pwp) <- c( ncol(coverages), ncol(coverages) )

for (i in 1:ncol(coverages)) {
    for (j in 1:ncol(coverages)) {
        i0 <- 2*i-1; i1 <- 2*i
        j0 <- 2*j-1; j1 <- 2*j
        weights <- sqrt(coverages[,i]) * sqrt(coverages[,j])
        if (i!=j) {
            probs <- counts[,i0]/coverages[,i] * counts[,j1]/coverages[,j] + counts[,i1]/coverages[,i] * counts[,j0]/coverages[,j]
            pwp[i,j] <- sum( weights*probs, na.rm=TRUE ) / sum(weights)
        } else {
            usethese <- ( coverages[,i] > 1 )
            probs <- 2 * counts[,i0] * counts[,i1] / ( coverages[,i] * (coverages[,i]-1) )
            pwp[i,i] <- sum( (weights*probs)[usethese] ) / sum(weights[usethese])
        }
    }
}

write.table(counts, col.names=FALSE, row.names=FALSE, sep=' ', file="test-pwp.counts")
awk <- scan(pipe("cat test-pwp.counts | awk -f ../weighted-pwp.awk"))

pwp.awk <- pwp
pwp.awk[upper.tri(pwp.awk,diag=TRUE)] <- awk
pwp.awk[lower.tri(pwp.awk)] <- t(pwp.awk)[lower.tri(pwp.awk)]

stopifnot( all( abs(pwp - pwp.awk) < 1e-6 ) )
