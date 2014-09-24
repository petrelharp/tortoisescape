#!/usr/bin/Rscript

source("resistance-fns.R")
require(raster)

layer.prefix <- commandArgs(TRUE)[1]

# for (layer.prefix in c( "../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_", "../geolayers/TIFF/10x/crop_resampled_masked_aggregated_10x_") ) {

    layer.names <- c("annual_precip","barren_30","bd_ss2_st_30","eastness_30","lat_gcs_30","lon_gcs_30")

    init.params <- c( beta=1.0, gamma=rep(.01,length(layer.names)), delta=rep(.01,length(layer.names)) )

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
    # generator matrix
    ij <- adjacent(onelayer,cells=nonmissing,target=nonmissing,directions=4,pairs=TRUE,sorted=TRUE) # to and from cells both loc
    ij <- ij[,2:1]
    stopifnot( all(ij[,1] != ij[,2]) ) ## NO DIAGONAL

    G <- sparseMatrix( i=ij[,1], j=ij[,2], x=1.0 )
    Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )

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
    # }
