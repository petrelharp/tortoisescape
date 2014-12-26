#!/usr/bin/Rscript

usage <- "
Load up everything already computed into one .RData file.  Usage:
    Rscript (layer prefix) (subdir) (layer file) (pairwise divergence file)
e.g.
    Rscript setup-inference.R ../geolayers/multigrid/512x/crm_ 512x six-raster-list ../pairwisePi/alleleCounts_1millionloci.pwp
where
    (layer prefix) = prefix to look for raster files in
    (subdir) = where to put stuff
    (layer file) = file with names of layers to use
    (pairwise divergence file) = file with UPPER triangle of matrix of pairwise divergences (including diagonals)
"

if (length(commandArgs(TRUE))<4) { stop(usage) }

if (!interactive()) {
    layer.prefix <- commandArgs(TRUE)[1]
    subdir <- commandArgs(TRUE)[2]
    layer.file <- commandArgs(TRUE)[3]
    pimat.file <- commandArgs(TRUE)[4]
} else {
    layer.prefix <- c("../geolayers/TIFF/500x/500x_")
    subdir <- "500x"
    layer.file <- "six-raster-list"
    pimat.file <- "../pairwisePi/alleleCounts_1millionloci.pwp"
}
layer.names <- scan(layer.file,what="char") 

source("resistance-fns.R")
require(raster)

outfile <- paste(subdir,"/",basename(layer.prefix),basename(layer.file),"-","setup.RData",sep='')

source.ls <- ls()

# get precomputed G
load(paste(subdir,"/",basename(layer.prefix),"_",basename(layer.file),"_G.RData",sep=''))
load(paste(subdir,"/",basename(layer.prefix),"_",basename(layer.file),"_nonmissing.RData",sep=''))
# Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )

###
# layer whatnot

layers <- do.call( cbind, lapply( layer.names, function (layer.name) {
        scale( values( raster(paste(layer.prefix,layer.name,sep='')) )[nonmissing] ) 
    } ) )
stopifnot(nrow(layers)==nrow(G))

# tortoise locations
load(paste(subdir,"/",basename(layer.prefix),basename(layer.file),"_tortlocs.RData",sep=''))
nind <- length(locs)
na.indiv <- which( is.na( locs ) )
locs <- locs[-na.indiv]

# and neighborhoods
load( paste( subdir, "/", basename(layer.prefix), basename(layer.file), "_neighborhoods.RData", sep='' ) ) # provides 'neighborhoods'
neighborhoods <- lapply(neighborhoods[-na.indiv],function (x) { x[!is.na(x)] })

# pairwise divergence values
pimat.vals <- scan(pimat.file) # has UPPER with diagonal
pimat <- numeric(nind^2)
dim(pimat) <- c(nind,nind)
pimat[upper.tri(pimat,diag=TRUE)] <- pimat.vals
pimat[lower.tri(pimat,diag=FALSE)] <- t(pimat)[lower.tri(pimat,diag=FALSE)]
pimat <- pimat[-na.indiv,-na.indiv]

# scale to actual pairwise divergence, and then by 1/mutation rate
pimat <- pimat * .018 * 1e8

save(list=setdiff(ls(),source.ls), file=outfile)

cat("Saved to ", outfile, " .\n")
