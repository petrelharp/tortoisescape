#!/usr/bin/Rscript

source("resistance-fns.R")
require(raster)
rasterOptions(tmpdir=".")

if (!interactive()) { 
    layer.prefix <- commandArgs(TRUE)[1] 
    layer.file <- commandArgs(TRUE)[2]
} else {
    layer.file <- "../inference/six-raster-list"
    layer.prefix <- c("../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_")
}

# for (layer.prefix in c( "../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_", "../geolayers/TIFF/10x/crop_resampled_masked_aggregated_10x_", "../geolayers/TIFF/masked/crop_resampled_masked_" ) ) {

    # # minimal list to make sure this works with update.G below
    # layer.names <- c("annual_precip","lon_gcs_30")
    layer.names <- paste( basename(layer.file), "-na", sep='' )  # this is TRUE in cells with missing data

    init.params <- c( beta=1.0, gamma=rep(.01,length(layer.names)), delta=rep(.01,length(layer.names)) )

    ###
    # layer whatnot

    # get info out of na layer
    nalayer <- raster(paste(dirname(layer.prefix),"/",layer.names[1],sep=''))
    n <- dim(nalayer)[2]; m <- dim(nalayer)[1]
    nmap <- matrix(1:n*m,nrow=n)
    all.locs <- cbind( i=as.vector(row(nmap)), j=as.vector(col(nmap)) )
    nonmissing <- which(values(nalayer)==0)

    ###
    # generator matrix
    ij <- adjacent(nalayer,cells=nonmissing,target=nonmissing,directions=4,pairs=TRUE,sorted=TRUE) # to and from cells both loc
    ij <- ij[,2:1]
    stopifnot( all(ij[,1] != ij[,2]) ) ## NO DIAGONAL

    G <- sparseMatrix( i=match(ij[,1],nonmissing), j=match(ij[,2],nonmissing), x=1.0 )
    Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )

    layers <- cbind( values(nalayer)[nonmissing] ) # should be all zero, but heck
    stopifnot(nrow(layers)==nrow(G))

    transfn <- exp
    valfn <- function (gamma) { ( rowSums( layers * gamma[col(layers)], na.rm=TRUE ) ) }

    ndelta <- ngamma <- length(layer.names)
    update.G <- function(params) {
        beta <- params[1]
        gamma <- params[1+(1:ngamma)]
        delta <- params[1+ngamma+(1:ndelta)]
        return( beta * transfn(valfn(gamma))[G@i+1L] * transfn( valfn(delta)[G@i+1L] + valfn(delta)[Gjj] ) )
    }

    G@x <- update.G(init.params)

    save( G, update.G, ndelta, ngamma, transfn, valfn, layers, file=paste(basename(layer.prefix),"G.RData",sep=''))
    save( nonmissing, file=paste(basename(layer.prefix),"nonmissing.RData",sep=''))
    # }
