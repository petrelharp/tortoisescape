nsites <- 100
counts <- cbind( x0=rpois(nsites,1), x1=rpois(nsites,1), y0=rpois(nsites,2), y1=rpois(nsites,2), z0=rpois(nsites,10), z1=rpois(nsites,10) )
coverages <- counts[,2*(1:(ncol(counts)/2))] + counts[,2*(1:(ncol(counts)/2))-1]

coverage.cutoff <- 12

covmat <- numeric( ncol(coverages)^2 )
dim(covmat) <- c( ncol(coverages), ncol(coverages) )
dimnames(covmat) <- list( c("x","y","z"), c("x","y","z") )

for (i in 1:ncol(coverages)) {
    for (j in 1:ncol(coverages)) {
        i0 <- 2*i-1; i1 <- 2*i
        j0 <- 2*j-1; j1 <- 2*j
        weights <- (coverages[,i] * ( coverages[,j] - if(i==j){1}else{0}))
        weights[ (coverages[,i]>coverage.cutoff) | (coverages[,j]>coverage.cutoff) ] <- 0
        freq1 <- counts[,i0]/coverages[,i]
        freq2 <- counts[,j0]/coverages[,j]
        covmat[i,j] <- sum( weights*freq1*freq2, na.rm=TRUE ) / sum(weights) - ( sum( weights*freq1, na.rm=TRUE ) / sum(weights) ) * ( sum( weights*freq2, na.rm=TRUE ) / sum(weights) )
    }
}

cat(c("x","y","z\n"),file="test-pwp.counts",sep=" ")
write.table(counts, col.names=FALSE, row.names=FALSE, sep=' ', file="test-pwp.counts", append=TRUE)
awk.df <- read.table(pipe("cat test-pwp.counts | awk -f ../weighted-covariance.awk"))

covmat.awk <- covmat
covmat.awk[] <- NA
covmat.awk[ cbind( match(awk.df[,1],rownames(covmat.awk)), match(awk.df[,2],rownames(covmat.awk)) ) ] <- awk.df[,3]
covmat.awk[is.na(covmat.awk)] <- t(covmat.awk)[is.na(covmat.awk)]

stopifnot( all( abs(covmat - covmat.awk) < 1e-6 ) )
