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

hasdata <- with(torts, UTM_Zone!="12N" )
ut <- upper.tri(pimat,diag=FALSE) 
usethese <- which( outer(hasdata,hasdata,"&") & ut )
require(colorspace)
tcols <- rainbow(240)

ibd.list <- lapply( 1:nind, function (k) lm( pimat[k,hasdata] ~ tort.dists[k,hasdata] ) )
ibd.coef <- sapply(ibd.list,coef)

var_col <- function (varname) { 1 + outer( torts[[varname]], torts[[varname]], "!=" ) }

torts$Adapter11 <- ( torts$Adapter == "BFIDT-011" )

for (varname in c("Sex","Cleaned","Tissue.Box","Adapter11","Lib.prepped.","Sent.to.Berk.")) {
    plot( tort.dists[usethese], pimat[usethese], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5, col=var_col(varname)[usethese], main=varname )
    if (is.null(locator(1))) { break }
}

et110 <- which( torts$EM_Tort_ID == "etort-110" )
et50 <- which( torts$EM_Tort_ID == "etort-50" )
other11 <- which( torts$EM_Tort_ID != "etort-110" & ! torts$EM_Tort_ID == "etort-50" & torts$Adapter11 ) 

plot( tort.dists[usethese], pimat[usethese], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5, col=var_col("Adapter11")[usethese], main=varname )
points( tort.dists[et110,], pimat[et110,], col='green' )
points( tort.dists[et50,], pimat[et50,], col='purple' )
points( tort.dists[other11,other11], pimat[other11,other11], col='blue' )
ad11 <- with(torts,outer(Adapter11,Adapter11,"|"))
xx <- tort.dists[ut & outer(hasdata,hasdata,"&") & ad11]
yy <- pimat[ut & outer(hasdata,hasdata,"&") & ad11]
abline( coef( lm( yy ~ xx ) ), col='red', lwd=2 )
xx <- tort.dists[ut & outer(hasdata,hasdata,"&") & !ad11]
yy <- pimat[ut & outer(hasdata,hasdata,"&") & !ad11]
abline( coef( lm( yy ~ xx ) ), lwd=2 )

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

