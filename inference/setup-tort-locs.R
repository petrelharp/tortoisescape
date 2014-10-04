#!/usr/bin/Rscript

source("resistance-fns.R")
require(raster)
rasterOptions(tmpdir=".")

require(parallel)
numcores<-as.numeric(scan(pipe("cat /proc/cpuinfo | grep processor | tail -n 1 | awk '{print $3}'")))+1

# for (layer.prefix in c( "../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_", "../geolayers/TIFF/10x/crop_resampled_masked_aggregated_10x_", "../geolayers/TIFF/masked/crop_resampled_masked_" ) ) {

if (!interactive()) { 
    layer.prefix <- commandArgs(TRUE)[1] 
    subdir <- commandArgs(TRUE)[2] 
    layer.file <- commandArgs(TRUE)[3] 
} else {
    layer.prefix <- "../geolayers/TIFF/500x/500x_"
    subdir <- "500x"
    layer.file <- "six-raster-list"
}

onelayer <- raster(paste(layer.prefix,"dem_30",sep=''))

###
# tortoise locations
load("../tort.coords.rasterGCS.Robj")  # provides tort.coords.rasterGCS
orig.locs <- cellFromXY( onelayer, tort.coords.rasterGCS )
# #  locs.ij <- k.to.ij(locs,n)  # this produces j's that are off by one (zero-based?)
# locs.ij <- cbind( colFromX(onelayer,tort.coords.rasterGCS), rowFromY(onelayer,tort.coords.rasterGCS) )
# stopifnot( all( orig.locs == cellFromRowCol(onelayer,locs.ij[,2],locs.ij[,1]) ) )


# ok, but we need indices in NONMISSING ones
load(paste(subdir,"/",basename(layer.prefix),"_",basename(layer.file),"_nonmissing.RData",sep=''))
locs <- match(orig.locs,nonmissing)

save( locs, file=paste(subdir,"/",basename(layer.prefix),"tortlocs.RData",sep='') )

###
# and distances from all nonmissing locations to torts
if (FALSE) {

    all.locs.dists <- sapply( 1:length(locs), function (k) {
                values( distanceFromPoints( onelayer, tort.coords.rasterGCS[k] ) )[nonmissing]
            } )

    save( all.locs.dists, file=paste(subdir,"/",basename(layer.prefix),"alllocs.RData",sep='') )
}

##
# Compute neighborhoods of each tortoise
#  to use in mean hitting time computations

ndist <- 15000  # 15 km

neighborhoods <- mclapply( seq_along(tort.coords.rasterGCS) , function (k) {
        d_tort <- distanceFromPoints( onelayer, tort.coords.rasterGCS[k] )
        match( Which( d_tort <= max(ndist,minValue(d_tort)), cells=TRUE, na.rm=TRUE ), nonmissing )
    }, mc.cores=numcores )

save(neighborhoods, file=paste( subdir, "/", basename(layer.prefix), "_", basename(layer.file), "_neighborhoods.RData", sep='' ) )

stopifnot( all( sapply( seq_along(locs), function(k) { locs[k] %in% neighborhoods[[k]] } ) ) )
