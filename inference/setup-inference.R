#!/usr/bin/Rscript

source("resistance-fns.R")
require(raster)

source.ls <- ls()

if (!interactive()) {
    layer.prefix <- commandArgs(TRUE)[1]
    subdir <- commandArgs(TRUE)[2]
    layer.file <- commandArgs(TRUE)[3]
} else {
    layer.prefix <- c("../geolayers/TIFF/500x/500x_")
    subdir <- "500x"
    layer.file <- "six-raster-list"
    # layer.names <- c("imperv_30", "agp_250", "m2_ann_precip", "avg_rough_30", "dem_30", "bdrock_ss2_st")
}
layer.names <- scan(layer.file,what="char") 

# get precomputed G
load(paste(subdir,"/",basename(layer.prefix),"G.RData",sep=''))
load(paste(subdir,"/",basename(layer.prefix),"nonmissing.RData",sep=''))
Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )

###
# layer whatnot

layers <- do.call( cbind, lapply( layer.names, function (layer.name) {
        scale( values( raster(paste(layer.prefix,layer.name,sep='')) )[nonmissing] ) 
    } ) )
stopifnot(nrow(layers)==nrow(G))

# tortoise locations
load(paste(subdir,"/",basename(layer.prefix),"tortlocs.RData",sep=''))
nind <- length(locs)
na.indiv <- which( is.na( locs ) )
locs <- locs[-na.indiv]

# pairwise divergence values
pimat.vals <- scan("../pairwisePi/alleleCounts_1millionloci.pwp") # has UPPER with diagonal
pimat <- numeric(nind^2)
dim(pimat) <- c(nind,nind)
pimat[upper.tri(pimat,diag=TRUE)] <- pimat.vals
pimat[lower.tri(pimat,diag=FALSE)] <- t(pimat)[lower.tri(pimat,diag=FALSE)]
pimat <- pimat[-na.indiv,-na.indiv]

# scale to actual pairwise divergence, and then by 1/mutation rate
pimat <- pimat * .018 * 1e8

save(list=setdiff(ls(),source.ls), file=paste(subdir,"/",basename(layer.file),"-",basename(layer.prefix),"setup.RData",sep='') )
