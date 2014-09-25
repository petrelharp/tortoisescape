#!/usr/bin/Rscript

source("resistance-fns.R")
require(raster)
rasterOptions(tmpdir=".")

# for (layer.prefix in c( "../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_", "../geolayers/TIFF/10x/crop_resampled_masked_aggregated_10x_", "../geolayers/TIFF/masked/crop_resampled_masked_" ) ) {

if (!interactive()) { layer.prefix <- commandArgs(TRUE)[1] }

onelayer <- raster(paste(layer.prefix,"annual_precip",sep=''))

###
# tortoise locations
load("../tort.coords.rasterGCS.Robj")  # provides tort.coords.rasterGCS
orig.locs <- cellFromXY( onelayer, tort.coords.rasterGCS )
# #  locs.ij <- k.to.ij(locs,n)  # this produces j's that are off by one (zero-based?)
# locs.ij <- cbind( colFromX(onelayer,tort.coords.rasterGCS), rowFromY(onelayer,tort.coords.rasterGCS) )
# stopifnot( all( orig.locs == cellFromRowCol(onelayer,locs.ij[,2],locs.ij[,1]) ) )


# ok, but we need indices in NONMISSING ones
load(paste(basename(layer.prefix),"nonmissing.RData",sep=''))
locs <- match(orig.locs,nonmissing)

save( locs, file=paste(basename(layer.prefix),"tortlocs.RData",sep='') )

###
# and distances from all nonmissing locations to torts
if (dim(onelayer)[1]<5000) {

    all.locs.dists <- sapply( 1:length(locs), function (k) {
                values( distanceFromPoints( onelayer, tort.coords.rasterGCS[k] ) )[nonmissing]
            } )

    save( all.locs.dists, file=paste(basename(layer.prefix),"alllocs.RData",sep='') )
}
