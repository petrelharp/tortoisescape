#!/usr/bin/Rscript

options(error=dump.frames)

if (length(commandArgs(TRUE)) < 2) {
    cat("Usage: \n Rscript counts-pca.R infile outfile \n\n")
}

infile <- commandArgs(TRUE)[1]
outfile <- commandArgs(TRUE)[2]
# infile <- "exampleOutput/alleleCounts100k.txt"
# outfile <- "exampleOutput/alleleCounts100k-covmat.txt"

fin <- if (infile=="-") { stdin() } else { file(infile,open="r") }
c1 <- scan(fin,nlines=1)
counts <- c(c1,scan(fin))
close(fin)
dim(counts) <- c(length(c1), length(counts)/length(c1))
# note this is TRANSPOSED

ninds <- nrow(counts)/2
majors <- counts[2*(1:ninds)-1,]
totals <- counts[2*(1:ninds)-1,] + counts[2*(1:ninds),]

w.covmat <- matrix( NA, nrow=ninds, ncol=ninds )
for (ii in 1:ninds) for (jj in ii:ninds) {
    fi <- (majors/totals)[ii,]
    fj <- (majors/totals)[jj,]
    w <- sqrt( totals[ii,] * totals[jj,] )
    w.covmat[jj,ii] <- w.covmat[ii,jj] <- weighted.mean( fi * fj, w=w ) - weighted.mean( fi, w=w ) * weighted.mean( fj, w=w )
}

fout <- if (outfile=="-") { stdout() } else { file(outfile,open="w") }
write.table( w.covmat, file=fout, quote=FALSE, row.names=FALSE, col.names=FALSE )
close(fout)

# also do another version
if (outfile != "-") {
    outfile.norm <- paste(outfile,".norm")

    n.covmat <- matrix( NA, nrow=ninds, ncol=ninds )
    p <- colSums(majors) / colSums(totals)
    for (ii in 1:ninds) for (jj in ii:ninds) {
        pi <- (majors/totals)[ii,] - p
        pj <- (majors/totals)[jj,] - p
        w <- sqrt( totals[ii,] * totals[jj,] )
        n.covmat[jj,ii] <- n.covmat[ii,jj] <- weighted.mean( pi * pj, w=w ) - weighted.mean( pi, w=w ) * weighted.mean( pj, w=w )
    }

    fout <- file(outfile.norm,open="w")
    write.table( n.covmat, file=fout, quote=FALSE, row.names=FALSE, col.names=FALSE )
    close(fout)
}


if (FALSE) {

    ## testing code
    ninds <- 4
    nsites <- 1e6
    coverage <- matrix( rpois( ninds * nsites, lambda=rep(1:ninds,each=nsites) ), nrow=nsites )
    # chol.covmat <- diag(ninds)
    chol.covmat <- chol( crossprod( matrix( rexp(ninds^2,rate=5), nrow=ninds ) ) )
    freqs <- 1/(1+exp( matrix( 10*rnorm(ninds*nsites), nrow=nsites ) %*% chol.covmat ) )
    freq.covmat <- cov(freqs)
    counts <- matrix( rbinom( ninds*nsites, size=coverage, prob=freqs ), nrow=nsites )
    counts.covmat <- cov(counts)

    adj.covmat <- counts.covmat / ( crossprod( coverage ) / nsites )

    freq.covmat - adj.covmat

    rbind( colMeans( coverage * freqs * (1-freqs) + coverage^2 * freqs^2 ) - colMeans( freqs * coverage )^2,
        diag(counts.covmat) )

    # this estimates offdiagonals of freq.covmat well, but underestimates diagonal
    w.covmat <- matrix( NA, nrow=ninds, ncol=ninds )
    for (ii in 1:ninds) for (jj in ii:ninds) {
        fi <- (counts/coverage)[,ii]
        fj <- (counts/coverage)[,jj]
        w <- sqrt( coverage[,ii] * coverage[,jj] )
        w.covmat[jj,ii] <- w.covmat[ii,jj] <- weighted.mean( fi * fj, w=w ) - weighted.mean( fi, w=w ) * weighted.mean( fj, w=w )
    }

}
