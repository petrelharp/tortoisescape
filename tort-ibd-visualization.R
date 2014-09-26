#!/usr/bin/Rscript

load("tort.coords.rasterGCS.Robj")

require(raster)
layer <- raster("geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_dem_30.gri")

torts <- read.csv("1st_180_torts.csv",header=TRUE)
nind <- nrow(torts)

tort.dist.table <- read.table("1st180_pairwise_distances_sorted_redundancy_removed.txt",header=TRUE)
tort.dists <- numeric(nind^2); dim(tort.dists) <- c(nind,nind)
tort.dists[ cbind( match(tort.dist.table$etort1,torts$EM_Tort_ID), match(tort.dist.table$etort2,torts$EM_Tort_ID) ) ] <- tort.dist.table$DISTANCE
tort.dists <- tort.dists + t(tort.dists)

pimat.vals <- scan("pairwisePi/alleleCounts_1millionloci.pwp") # has UPPER with diagonal
pimat <- numeric(nind^2)
dim(pimat) <- c(nind,nind)
pimat[upper.tri(pimat,diag=TRUE)] <- pimat.vals
pimat[lower.tri(pimat,diag=FALSE)] <- t(pimat)[lower.tri(pimat,diag=FALSE)]

pcs <- read.csv("covmat/tort-PCs.csv")
pc.cols <- ifelse( pcs$PC1 > 0, "blue", "purple" )

for (k in 1:nind) { 
    png( file=paste("pngs/",torts$EM_Tort_ID[k],"-ibd.png",sep=''), width=8*144, height=4*144, pointsize=10, res=144 )
    layout(t(1:2))
    plot(layer)
    points(tort.coords.rasterGCS,pch=20,cex=1,col=pc.cols)
    points(tort.coords.rasterGCS[k],cex=2,col='red')
    plot( tort.dists[upper.tri(pimat)], pimat[upper.tri(pimat)], pch=20, cex=.5, col=adjustcolor("black",.25) )
    points( tort.dists[k,], pimat[k,], pch=20, col=pc.cols )
    # if (is.null(locator(1))) { break }
    dev.off()
}
