#!/usr/bin/Rscript

angsd.pimat.vals <- scan("350000.pwp")  # has LOWER triangle of entries, without diagonal
robust.pimat.vals <- scan("alleleCounts_1millionloci.pwp") # robust version, has UPPER with diagonal
nind <- 180

robust.pimat <- angsd.pimat <- numeric(nind^2)
dim(robust.pimat) <- dim(angsd.pimat) <- c(nind,nind)
angsd.pimat[lower.tri(angsd.pimat,diag=FALSE)] <- angsd.pimat.vals
angsd.pimat <- angsd.pimat + t(angsd.pimat)
robust.pimat[upper.tri(pimat,diag=TRUE)] <- robust.pimat.vals
robust.pimat[lower.tri(pimat,diag=FALSE)] <- t(robust.pimat)[lower.tri(pimat,diag=FALSE)]

mean(angsd.pimat)
mean(robust.pimat)

torts <- read.csv("../1st_180_torts.csv",header=TRUE)

# pairwise distances
tort.dist.table <- read.table("../1st180_pairwise_distances_sorted_redundancy_removed.txt",header=TRUE)
tort.dists <- numeric(nind^2); dim(tort.dists) <- c(nind,nind)
tort.dists[ cbind( match(tort.dist.table$etort1,torts$EM_Tort_ID), match(tort.dist.table$etort2,torts$EM_Tort_ID) ) ] <- tort.dist.table$DISTANCE
tort.dists <- tort.dists + t(tort.dists)

# coverage
coverages <- read.csv("../coverage_info.csv")
torts$coverage <- coverages$sequence_yield_megabases

###
# compare two methods
for (meth in c("angsd","robust")) {
    pimat <- get(paste(meth,"pimat",sep='.'))

    ut <- upper.tri(pimat,diag=FALSE) 
    require(colorspace)
    tcols <- rainbow(240)
    lcols <- adjustcolor(rainbow(nlevels(torts$Location_ID)),0.75)

    ibd.list <- lapply( 1:nind, function (k) lm( pimat[k,] ~ tort.dists[k,] ) )
    ibd.coef <- sapply(ibd.list,coef)

    var_col <- function (varname) { 1 + outer( torts[[varname]], torts[[varname]], "!=" ) }

    torts$Adapter11 <- ( torts$Adapter == "BFIDT-011" )

    if (FALSE) {
        for (varname in c("Sex","Cleaned","Tissue.Box","Adapter11","Lib.prepped.","Sent.to.Berk.")) {
            plot( tort.dists[ut], pimat[ut], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5, col=var_col(varname)[ut], main=varname )
            if (is.null(locator(1))) { break }
        }
    }

    et110 <- which( torts$EM_Tort_ID == "etort-110" )
    et50 <- which( torts$EM_Tort_ID == "etort-50" )
    other11 <- which( torts$EM_Tort_ID != "etort-110" & ! torts$EM_Tort_ID == "etort-50" & torts$Adapter11 ) 

    pdf(file=paste(meth,"pairwise-pi.pdf",sep='-'), width=10, height=6, pointsize=10)

    plot( tort.dists[ut], pimat[ut], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5 )
    abline( coef( lm( pimat[ut] ~ tort.dists[ut] ) ) ) 

    plot( tort.dists[ut], pimat[ut], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5, col=lcols[torts$Location_ID][row(pimat)[ut]] ) # tcols[row(pimat)[ut]] )
    abline( coef( lm( pimat[ut] ~ tort.dists[ut] ) ) ) 
    legend("bottomright",pch=20,col=lcols,legend=levels(torts$Location_ID),cex=0.5)

    plot( tort.dists[ut], pimat[ut], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5, col=lcols[torts$Location_ID][col(pimat)[ut]] ) # tcols[col(pimat)[ut]] )
    abline( coef( lm( pimat[ut] ~ tort.dists[ut] ) ) ) 
    legend("bottomright",pch=20,col=lcols,legend=levels(torts$Location_ID),cex=0.5)

    layout(matrix(1:4,nrow=2))
    plot( torts$coverage[row(pimat)[ut]], pimat[ut], cex=0.5, pch=20, main='coverage 1' )
    plot( torts$coverage[col(pimat)[ut]], pimat[ut], cex=0.5, pch=20, main='coverage 2' )
    plot( torts$coverage[col(pimat)[ut]]+torts$coverage[row(pimat)[ut]], pimat[ut], cex=0.5, pch=20, main='mean coverage' )
    plot( torts$coverage[col(pimat)[ut]]-torts$coverage[row(pimat)[ut]], pimat[ut], cex=0.5, pch=20, main='diff coverage' )

    layout(1)
    plot( tort.dists[ut], pimat[ut], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5, col=var_col("Adapter11")[ut] )
    points( tort.dists[et110,], pimat[et110,], col='green' )
    points( tort.dists[et50,], pimat[et50,], col='purple' )
    points( tort.dists[other11,other11], pimat[other11,other11], col='blue' )
    ad11 <- with(torts,outer(Adapter11,Adapter11,"|"))
    xx <- tort.dists[ut & ad11]
    yy <- pimat[ut & ad11]
    abline( coef( lm( yy ~ xx ) ), col='red', lwd=2 )
    xx <- tort.dists[ut & !ad11]
    yy <- pimat[ut & !ad11]
    abline( coef( lm( yy ~ xx ) ), lwd=2 )
    legend("bottomright", pch=c(20,1,1,1), col=c("red","green","purple","blue"), legend=c("adapter 11-others", "etort-110", "etort-50", "adapter 11-adapter 11") )

    dev.off()

}
