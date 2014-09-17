#!/usr/bin/Rscript

pimat.vals <- scan("100000.cov")  # has upper triangle of entries
nind <- 180

pimat <- numeric(nind^2)
dim(pimat) <- c(nind,nind)
pimat[lower.tri(pimat,diag=FALSE)] <- pimat.vals
pimat <- pimat + t(pimat)

mean(pimat)

torts <- read.csv("../1st_180_torts.csv",header=TRUE)

tort.dists <- sqrt( with(torts, outer(Easting,Easting,"-")^2 + outer(Northing,Northing,"-")^2 ) )

hasdata <- with(torts, UTM_Zone=="11N" )
usethese <- which( outer(hasdata,hasdata,"&") & upper.tri(pimat,diag=FALSE) )
require(colorspace)
tcols <- rainbow(240)

ibd.list <- lapply( 1:nind, function (k) lm( pimat[k,hasdata] ~ tort.dists[k,hasdata] ) )
ibd.coef <- sapply(ibd.list,coef)

pdf(file="pairwise-pi-first-100000.pdf",width=10, height=6, pointsize=10)

plot( tort.dists[usethese], pimat[usethese], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5 )
abline( coef( lm( pimat[usethese] ~ tort.dists[usethese] ) ) ) 

plot( tort.dists[usethese], pimat[usethese], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5, col=tcols[row(pimat)[usethese]] )
abline( coef( lm( pimat[usethese] ~ tort.dists[usethese] ) ) ) 

plot( tort.dists[usethese], pimat[usethese], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5, col=tcols[col(pimat)[usethese]] )
abline( coef( lm( pimat[usethese] ~ tort.dists[usethese] ) ) ) 

plot( tort.dists[usethese], pimat[usethese], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5, col=tcols[col(pimat)[usethese]] )
for (k in (1:nind)[hasdata]) {
    abline( coef( lm( pimat[k,hasdata] ~ tort.dists[k,hasdata] ) ), col=tcols[k], lwd=2 )
}

dev.off()

