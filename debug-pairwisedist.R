#!/usr/bin/Rscript
require(maps)
require(maptools)
require(raster)


###
# version 1
torts <- read.csv("1st_180_torts.csv",header=TRUE)
nind <- nrow(torts)

# pairwise distances
tort.dist.table <- read.table("1st180_pairwise_distances_sorted_redundancy_removed.txt",header=TRUE)
tort.dists <- numeric(nind^2); dim(tort.dists) <- c(nind,nind)
tort.dists[ cbind( match(tort.dist.table$etort1,torts$EM_Tort_ID), match(tort.dist.table$etort2,torts$EM_Tort_ID) ) ] <- tort.dist.table$DISTANCE
tort.dists <- tort.dists + t(tort.dists)


###
# version 2
layer.prefix <- "geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_"
onelayer <- raster(paste(layer.prefix,"annual_precip",sep=''))
nonmissing <- which( !is.na(values(onelayer)) )

load("tort.coords.rasterGCS.Robj")  # provides tort.coords.rasterGCS
orig.locs <- cellFromXY( onelayer, tort.coords.rasterGCS )
locs <- match(orig.locs,nonmissing)
all.locs.dists <- sapply( 1:length(locs), function (k) {
            values( distanceFromPoints( onelayer, tort.coords.rasterGCS[k] ) )[nonmissing]
        } )

###
# compare

plot(onelayer)
points(tort.coords.rasterGCS)

plot( tort.dists[upper.tri(tort.dists)], all.locs.dists[locs,][upper.tri(tort.dists)] )
