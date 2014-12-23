#!/usr/bin/Rscript

usage <- "
Make a layer which records which locations are present (not missing) in common across the layers given.  Usage:
    Rscript make-overlap-na-layer.R (layer prefix) (layer file)
e.g.
    Rscript make-overlap-na-layer.R ../geolayers/multigrid/512x/crm_ six-raster-list
where
    (layer prefix) = prefix to look for raster files in
    (layer file) = file with names of layers to use

"

if (length(commandArgs(TRUE))<2) { stop(usage) }

if (!interactive()) {
    layer.prefix <- commandArgs(TRUE)[1]
    layer.file <- commandArgs(TRUE)[2]
} else {
    # layer.prefix <- c("TIFF/masked/crop_resampled_masked_")
    layer.prefix <- c("../geolayers/TIFF/500x/500x_")
    layer.file <- "../inference/twentyfour-raster-list"
}

require(raster)
rasterOptions(tmpdir=".")

layer.files <- list.files(dirname(layer.prefix),paste(basename(layer.prefix),".*gri",sep=''),full.names=TRUE)
layer.names <- gsub(layer.prefix,"",gsub(".gri","",layer.files,fixed=TRUE))
use.ref <- scan(layer.file,what="char")
use.ref <- intersect( use.ref, layer.names )

use.files <- layer.files[match(use.ref,layer.names)]
use.names <- layer.names[match(use.ref,layer.names)]

names(use.files) <- use.names

refname <- grep("dem_30",layer.names)[1]
reflayer <- is.na( raster(layer.files[refname]) ) # dem

for (lf in use.files) {
    other <- raster(lf) 
    reflayer <- reflayer | is.na(other)
}

outfile <- paste(dirname(layer.prefix),"/",basename(layer.file),"-na",sep='')
cat('writing to: ', outfile,' .\n')
writeRaster( reflayer, file=outfile, overwrite=TRUE )

if (FALSE) {
    require(parallel)
    refname <- grep("annual_precip",layer.names)
    reflayer <- raster(layer.files[refname]) #annual_precip

    subset.na <- mclapply( use.files, function (lf) {
            other <- raster(lf) 
            table( values( is.na(reflayer) ), values( is.na(other) ) )
        }, mc.cores=16 )

}
