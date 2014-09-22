#!/usr/bin/Rscript

source("resistance-fns.R")
require(raster)

layer.prefix <- c("../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_")
layer.names <- c("annual_precip","barren_30","bd_ss2_st_30","eastness_30","lat_gcs_30","lon_gcs_30")

init.params <- c( rep(.01,length(layer.names)), rep(.01,length(layer.names)) )

###
# layer whatnot

# get info out of one layer
onelayer <- raster(paste(layer.prefix,layer.names[1],sep=''))
n <- dim(onelayer)[2]; m <- dim(onelayer)[1]
nmap <- matrix(1:n*m,nrow=n)
all.locs <- cbind( i=as.vector(row(nmap)), j=as.vector(col(nmap)) )
onevals <- values(onelayer)
nonmissing <- which(!is.na(onevals))

layers <- sapply(layer.names, function (ll) {
            rast <- raster(paste(layer.prefix,ll,sep=''))
            # note this is ROW-ORDERED
            # so to plot do:  dim(x) <- dim(rast)[2:1]; image(x)
            vrast <- scale( values(rast)[nonmissing] )
            return(vrast)
        } )

###
# tortoise locations
load("../tort.coords.rasterGCS.Robj")  # provides tort.coords.rasterGCS
locs <- cellFromXY( onelayer, tort.coords.rasterGCS )
#  locs.ij <- k.to.ij(locs,n)  # this produces j's that are off by one (zero-based?)
locs.ij <- cbind( colFromX(onelayer,tort.coords.rasterGCS), rowFromY(onelayer,tort.coords.rasterGCS) )
stopifnot( all( locs == cellFromRowCol(onelayer,locs.ij[,2],locs.ij[,1]) ) )


###
# generator matrix
ij <- adjacent(onelayer,cells=nonmissing,target=nonmissing,directions=4,pairs=TRUE) # to and from cells both loc

transfn <- exp
valfn <- function (gamma) { ( rowSums( layers * gamma[col(layers)], na.rm=TRUE ) ) }

update.G <- function(params) {
    gamma <- params[1:length(delta)]
    delta <- params[length(delta)+(1:length(gamma))]
    G <- grid.adjacency(n,m,diag=FALSE,symmetric=FALSE)
    dp <- diff(G@p)
    jj <- rep(seq_along(dp),dp)
    G@x <- transfn(valfn(gamma))[G@i+1L] * transfn( valfn(delta)[G@i+1L] + valfn(delta)[jj] )
    diag(G) <- (-1) * rowSums(G)
    return(G)
}

