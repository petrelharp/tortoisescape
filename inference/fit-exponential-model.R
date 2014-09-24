#!/usr/bin/Rscript

source("resistance-fns.R")
require(raster)

layer.prefix <- c("../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_")

# get precomputed G
load(paste(basename(layer.prefix),"G.RData",sep=''))
Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )

###
# layer whatnot

init.params <- c( beta=1.0, gamma=rep(.01,length(layer.names)), delta=rep(.01,length(layer.names)) )

layer.names <- c("annual_precip","barren_30","eastness_30","lat_gcs_30","lon_gcs_30")
layers <- sapply(layer.names, function (ll) {
            rast <- raster(paste(layer.prefix,ll,sep=''))
            # note this is ROW-ORDERED
            # so to plot do:  dim(x) <- dim(rast)[2:1]; image(x)
            vrast <- scale( values(rast)[nonmissing] )
            return(vrast)
        } )
stopifnot(nrow(layers)==nrow(G))

G@x <- update.G(init.params)

###
# tortoise locations
load("../tort.coords.rasterGCS.Robj")  # provides tort.coords.rasterGCS
locs <- cellFromXY( onelayer, tort.coords.rasterGCS )
#  locs.ij <- k.to.ij(locs,n)  # this produces j's that are off by one (zero-based?)
locs.ij <- cbind( colFromX(onelayer,tort.coords.rasterGCS), rowFromY(onelayer,tort.coords.rasterGCS) )
stopifnot( all( locs == cellFromRowCol(onelayer,locs.ij[,2],locs.ij[,1]) ) )

# get some initial values for the jacobi iterative solver

init.hts <- matrix(1,nrow=nrow(G),ncol=length(locs))
system.time( jacobi.true.hts <- hitting.jacobi(locs,G,init.hts) )


