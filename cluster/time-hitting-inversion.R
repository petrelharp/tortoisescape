#!/usr/bin/Rscript

source("resistance-fns.R")
require(raster)

layer.prefix <- c("../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_")

load(paste(basename(layer.prefix),"G.RData",sep=''))


# use function "adjacent" (returns )
system.time(ij <- adjacent(rast,cells=loc,target=loc)) # to and from cells both loc

# clear rast and loc from memory
rm(rast,loc)
gc()

# declare matrix
system.time(G <- sparseMatrix(ij[,1],ij[,2],x=rnorm(length(ij[,1]),0,1)))
