#!/usr/bin/Rscript

# get weighted distances
# between pairs of torts
# summing across environmental layers
# along intervening raster pixels


source("resistance-fns.R")
require(raster)
rasterOptions(tmpdir=".")

require(parallel)
numcores<-as.numeric(scan(pipe("cat /proc/cpuinfo | grep processor | tail -n 1 | awk '{print $3}'")))+1

# for (layer.prefix in c( "../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_", "../geolayers/TIFF/10x/crop_resampled_masked_aggregated_10x_", "../geolayers/TIFF/masked/crop_resampled_masked_" ) ) {

if (!interactive()) { 
    layer.prefix <- commandArgs(TRUE)[1] 
    subdir <- commandArgs(TRUE)[2] 
    layer.name <- commandArgs(TRUE)[3]
} else {
    layer.prefix <- "../geolayers/TIFF/500x/500x_"
    subdir <- "500x"
    layer.name <- "annual_precip"
}

layer.file <- paste(layer.prefix,layer.name,sep='')
layer <- raster(layer.file)

###
# tortoise locations
load("../tort.coords.rasterGCS.Robj")  # provides tort.coords.rasterGCS
orig.locs <- cellFromXY( layer, tort.coords.rasterGCS )
tort.coords.raw <- coordinates(tort.coords.rasterGCS)

### tort info
torts <- read.csv("../1st_180_torts.csv",header=TRUE)
nind <- nrow(torts)

# pairwise divergence values
pimat.vals <- scan("../pairwisePi/alleleCounts_1millionloci.pwp") # has UPPER with diagonal
pimat <- numeric(nind^2)
dim(pimat) <- c(nind,nind)
pimat[upper.tri(pimat,diag=TRUE)] <- pimat.vals
pimat[lower.tri(pimat,diag=FALSE)] <- t(pimat)[lower.tri(pimat,diag=FALSE)]

# pairwise distances
tort.dist.table <- read.table("../1st180_pairwise_distances_sorted_redundancy_removed.txt",header=TRUE)
tort.dists <- numeric(nind^2); dim(tort.dists) <- c(nind,nind)
tort.dists[ cbind( match(tort.dist.table$etort1,torts$EM_Tort_ID), match(tort.dist.table$etort2,torts$EM_Tort_ID) ) ] <- tort.dist.table$DISTANCE
tort.dists <- tort.dists + t(tort.dists)


# pairwise paths between tortoises
tort.pairs <- combn( 1:nrow(torts), 2 )

env.dists <- unlist( mclapply( 1:ncol(tort.pairs), function (k) {
            t1 <- tort.pairs[1,k]
            t2 <- tort.pairs[2,k]
            x <- cellFromLine( layer, SpatialLines( list( Lines( list(Line( tort.coords.raw[c(t1,t2),] )), ID=paste(torts$EM_Tort_ID[c(t1,t2)],collapse='-') ) ), 
                            proj4string=crs(proj4string(tort.coords.rasterGCS)) ) )
            sum( abs(values(layer)[x[[1]]]) )
        }, mc.cores=numcores ) )

edists <- data.frame( tort1=torts$EM_Tort_ID[tort.pairs[1,]], tort2=torts$EM_Tort_ID[tort.pairs[2,]], edist=env.dists )

write.table( edists, file=paste(subdir,"/", layer.name, "-edist.tsv", sep=''), sep='\t', row.names=FALSE, quote=FALSE )

# NOT WORKING:
# tort.cells <- mclapply( cellFromLine( layer, SpatialLines( lapply( 1:ncol(tort.pairs), function (k) {
#                 t1 <- tort.pairs[1,k]
#                 t2 <- tort.pairs[2,k]
#                 Lines( list(Line( tort.coords.raw[c(t1,t2),] )), ID=paste(torts$EM_Tort_ID[c(t1,t2)],collapse='-') )
#             }), proj4string=crs(proj4string(tort.coords.rasterGCS)) ) ), 
#         function (x) sum( abs(values(layer)[x]) ), 
#     , mc.cores=numcores  )

if (FALSE) {
    tort.SpatialLines <- SpatialLines( lapply( 1:ncol(tort.pairs), function (k) {
                 t1 <- tort.pairs[1,k]
                 t2 <- tort.pairs[2,k]
                 Lines( list(Line( tort.coords.raw[c(t1,t2),] )), ID=paste(torts$EM_Tort_ID[c(t1,t2)],collapse='-') )
             }), proj4string=crs(proj4string(tort.coords.rasterGCS)) )
    plot(layer)
    lines(tort.SpatialLines)
}
