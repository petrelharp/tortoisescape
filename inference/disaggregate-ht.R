#!/usr/bin/Rscript

##
# create matrix to convert hitting times on one grid to another


source("resistance-fns.R")
require(raster)

if (!interactive()) {
    layer.prefix.1 <- commandArgs(TRUE)[1]
    layer.prefix.2 <- commandArgs(TRUE)[2]
    prev.ht <- commandArgs(TRUE)[3]
} else {
    layer.prefix.1 <- "../geolayers/TIFF/500x/500x_"
    layer.prefix.2 <- "../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_"
    subdir <- "500x"
    prev.ht <- "500x/six-raster-list-hitting-times.tsv"
}

layer.name <- "dem_30"
layer.1 <- raster(paste(layer.prefix.1,layer.name,sep=''))
layer.2 <- raster(paste(layer.prefix.2,layer.name,sep=''))


