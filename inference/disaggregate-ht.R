#!/usr/bin/Rscript

##
# convert hitting times on one grid to another


source("resistance-fns.R")
require(raster)

require(parallel)
numcores<-as.numeric(scan(pipe("cat /proc/cpuinfo | grep processor | tail -n 1 | awk '{print $3}'")))+1


if (!interactive()) {
    layer.prefix.1 <- commandArgs(TRUE)[1]
    layer.prefix.2 <- commandArgs(TRUE)[2]
    subdir.1 <- commandArgs(TRUE)[3]
    subdir.2 <- commandArgs(TRUE)[4]
    prev.ht <- commandArgs(TRUE)[5]
} else {
    layer.prefix.1 <- "../geolayers/TIFF/500x/500x_"
    layer.prefix.2 <- "../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_"
    subdir.1 <- "500x"
    subdir.2 <- "100x"
    prev.ht <- "500x/six-raster-list-hitting-times.tsv"
}

hts <- read.table(prev.ht,header=TRUE)

env <- new.env()
load( paste(subdir.1, "/", basename(layer.prefix.1),"nonmissing.RData",sep=''), envir=env ) # provides nonmissing
nonmissing.1 <- with( env, nonmissing )
load( paste(subdir.2, "/", basename(layer.prefix.2),"nonmissing.RData",sep=''), envir=env ) # provides nonmissing
nonmissing.2 <- with( env, nonmissing )

layer.name <- "dem_30"
layer.1 <- raster(paste(layer.prefix.1,layer.name,sep=''))
layer.2 <- raster(paste(layer.prefix.2,layer.name,sep=''))

values(layer.1)[-nonmissing.1] <- NA
values(layer.2)[-nonmissing.2] <- NA

##

checkit <- TRUE

new.hts <- do.call( rbind, mclapply( 1:ncol(hts), function (k) {
        values(layer.1)[nonmissing.1] <- hts[,k]
        layer.1.dis <- crop( disaggregate( layer.1, fact=5, method='bilinear' ), layer.2 )
        stopifnot( all( dim(layer.1.dis)==dim(layer.2) ) )
        # can skip this step, hopefully
        if (checkit) {
            layer.1.dis.res <- resample( layer.1.dis, layer.2 )
            stopifnot( all( abs( values(layer.1.dis)[nonmissing.2] - values(layer.1.dis.res)[nonmissing.2] ) < 1e-8 ) )
        }
        # get values out
        return( values(layer.1.dis)[nonmissing.2] )
    }, mc.cores=numcores ) )

if (FALSE) {
    layout(t(1:4))
    plot(layer.1)
    plot(layer.1.dis)
    plot(layer.1.dis.res)
    plot(layer.1.dis.res-layer.1.dis)
    # plot(layer.2)
}
