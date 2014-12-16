#!/usr/bin/Rscript

usage <- "
Find the cell in which each tortoise falls, and precompute their neighborhood.  Usage:
    Rscript (layer prefix) (subdir) (layer file)
e.g.
    Rscript setup-tort-locs.R ../geolayers/multigrid/512x/crm_ 512x six-raster-list
where
    (layer prefix) = prefix to look for raster files in
    (subdir) = where to put results
    (layer file) = file with names of layers to use

"

if (length(commandArgs(TRUE))<3) { stop(usage) }

if (!interactive()) { 
    layer.prefix <- commandArgs(TRUE)[1] 
    subdir <- commandArgs(TRUE)[2] 
    layer.file <- commandArgs(TRUE)[3] 
} else {
    layer.prefix <- "../geolayers/TIFF/500x/500x_"
    subdir <- "500x"
    layer.file <- "six-raster-list"
}

source("resistance-fns.R")
require(raster)
rasterOptions(tmpdir=".")

require(parallel)
numcores<-getcores()

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

locs.outfile <- paste(subdir,"/",basename(layer.prefix),"tortlocs.RData",sep='') 
save( locs, file=locs.outfile )
cat("Saved to ", locs.outfile, " .\n")

if (FALSE) {
    ###
    # and distances from all nonmissing locations to torts
    all.locs.dists <- sapply( 1:length(locs), function (k) {
                values( distanceFromPoints( onelayer, tort.coords.rasterGCS[k] ) )[nonmissing]
            } )

    save( all.locs.dists, file=paste(subdir,"/",basename(layer.prefix),"alllocs.RData",sep='') )
}

##
# Compute neighborhoods of each tortoise
#  to use in mean hitting time computations

ndist <- 15000  # 15 km

neighborhoods <- get.neighborhoods( ndist=ndist, locations=tort.coords.rasterGCS, nonmissing=nonmissing, layer=onelayer, numcores=numcores, na.rm=FALSE )
# mclapply( seq_along(tort.coords.rasterGCS) , function (k) {
#         d_tort <- distanceFromPoints( onelayer, tort.coords.rasterGCS[k] )
#         match( Which( d_tort <= max(ndist,minValue(d_tort)), cells=TRUE, na.rm=TRUE ), nonmissing )
#     }, mc.cores=numcores )

outfile <- paste( subdir, "/", basename(layer.prefix), basename(layer.file), "_neighborhoods.RData", sep='' )
save(neighborhoods, file=outfile )

cat("Saved to ", outfile, " .\n")

stopifnot( all( sapply( seq_along(locs), function(k) { locs[k] %in% neighborhoods[[k]] } ) ) )


