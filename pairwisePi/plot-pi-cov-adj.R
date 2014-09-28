#!/usr/bin/Rscript
#require(rgdal)
require(maps)
require(maptools)
require(raster)

angsd.pimat.vals <- scan("350000.pwp")  # has LOWER triangle of entries, without diagonal
robust.pimat.vals <- scan("alleleCounts_1millionloci.pwp") # robust version, has UPPER with diagonal
nind <- 180

robust.pimat <- angsd.pimat <- numeric(nind^2)
dim(robust.pimat) <- dim(angsd.pimat) <- c(nind,nind)
angsd.pimat[lower.tri(angsd.pimat,diag=FALSE)] <- angsd.pimat.vals
angsd.pimat <- angsd.pimat + t(angsd.pimat)
robust.pimat[upper.tri(robust.pimat,diag=TRUE)] <- robust.pimat.vals
robust.pimat[lower.tri(robust.pimat,diag=FALSE)] <- t(robust.pimat)[lower.tri(robust.pimat,diag=FALSE)]

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

# elevation raster
#elev.file <- "../geolayers/TIFF/10x/crop_resampled_masked_aggregated_10x_dem_30.gri"
#elev <- raster(elev.file)
# and tortoise locs on it
#load("../tort.coords.rasterGCS.Robj")
#load("../county_lines.Robj")  # @gbradburd: how was this produced?


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

    # subtract off effect of coverage on Pi, then plot (actually both at the same time)

    pdf(file=paste(meth,"pairwise-pi-cov-adj.pdf",sep='-'), width=10, height=6, pointsize=10)

	# find trendline in pi vs. coverage
	cov.adj <- coef( lm( pimat[ut] ~ torts$coverage[col(pimat)[ut]]+torts$coverage[row(pimat)[ut]] ) )

	# find mean
	mean_cov <- mean(torts$coverage[col(pimat)[ut]]+torts$coverage[row(pimat)[ut]])

#	# adjust values
	pimat.adj <- pimat[ut] - cov.adj[[2]]*(torts$coverage[col(pimat)[ut]]+torts$coverage[row(pimat)[ut]]-mean_cov) # works fine, just don't subset it twice

	# plot adjusted values in red
	plot( tort.dists[ut], pimat.adj, xlab="geographic distance", ylab="coverage adjusted pairwise divergence", pch=20, cex=0.5, col=adjustcolor("red",.75) )
	abline( coef( lm( pimat.adj ~ tort.dists[ut] ) ), col=adjustcolor("red",.75) ) # this seems to work

	#plot original values in black
	matplot( tort.dists[ut], pimat[ut], pch=20, cex=0.5, col=adjustcolor("black",.75), add=T )
	abline( coef( lm( pimat[ut] ~ tort.dists[ut] ) ),col=adjustcolor("black",.75) ) 

	#plot lines between old and new
	matplot( t(cbind(tort.dists[ut],tort.dists[ut])), t(cbind(pimat[ut],pimat.adj)),type="l",col=adjustcolor("black",.2), lty=1, add=T )

    dev.off()

}
