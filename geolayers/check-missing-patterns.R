#!/usr/bin/Rscript

if (!interactive()) {
    layer.prefix <- commandArgs(TRUE)[1]
    layer.file <- commandArgs(TRUE)[2]
} else {
    layer.file <- "../inference/twentyfour-raster-list"
    # layer.prefix <- c("TIFF/masked/crop_resampled_masked_")
    layer.prefix <- c("TIFF/100x/crop_resampled_masked_aggregated_100x_")
}

require(raster)


layer.files <- list.files(dirname(layer.prefix),paste(basename(layer.prefix),".*gri",sep=''),full.names=TRUE)
layer.names <- gsub(layer.prefix,"",gsub(".gri","",layer.files,fixed=TRUE))
use.ref <- scan(layer.file,what="char")
use.ref <- intersect( use.ref, layer.names )

use.files <- layer.files[match(use.ref,layer.names)]
use.names <- layer.names[match(use.ref,layer.names)]

names(use.files) <- use.names

refname <- grep("annual_precip",layer.names)
reflayer <- is.na( raster(layer.files[refname]) ) #annual_precip

for (lf in use.files) {
    other <- raster(lf) 
    reflayer <- reflayer | is.na(other)
}

writeRaster( reflayer, file=paste(dirname(layer.prefix),"/",basename(layer.file),"-na",sep='') )

if (FALSE) {
    require(parallel)
    refname <- grep("annual_precip",layer.names)
    reflayer <- raster(layer.files[refname]) #annual_precip

    subset.na <- mclapply( use.files, function (lf) {
            other <- raster(lf) 
            table( values( is.na(reflayer) ), values( is.na(other) ) )
        }, mc.cores=16 )

}
