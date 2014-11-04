#!/usr/bin/Rscript

usage <- '
    Get hitting times with a list of landscape layers:
        Rscript multiscale-hitting-times.R (layer prefix) (subdir) (layer file) (parameter file) (method) [initial guess] [max running time] [output file]
      e.g.
        Rscript multiscale-hitting-times.R ../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_ 100x ../inference/six-raster-list simple-init-params-six-raster-list.tsv CG

    Here `method` is either "analytic" or "CG".
'

if (!interactive()) {
    if (length(commandArgs(TRUE))<5) { cat(usage) }
    layer.prefix <- commandArgs(TRUE)[1]
    subdir <- commandArgs(TRUE)[2]
    layer.file <- commandArgs(TRUE)[3]
    param.file <- commandArgs(TRUE)[4] 
    method <- commandArgs(TRUE)[5] 
    prev.ht <- if (length(commandArgs(TRUE))>5) { commandArgs(TRUE)[6] } else { NULL } 
    maxtime <- if (length(commandArgs(TRUE))>6) { commandArgs(TRUE)[7] } else { 6*60*60 } 
    outfile <- if (length(commandArgs(TRUE))>7) { commandArgs(TRUE)[8] } else { NULL }
} else {
    layer.prefix <- "../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_"
    subdir <- "100x"
    layer.file <- "../inference/six-raster-list"
    param.file <- "simple-init-params-six-raster-list.tsv"
    method <- "CG"
    prev.ht <- "100x/500x-aggregated-hitting-times.tsv"
    maxtime <- 2*60
    outfile <- NULL

    # layer.prefix <- "../geolayers/TIFF/500x/500x_"
    # subdir <- "500x"
    # layer.file <- "../inference/six-raster-list"
    # param.file <- "simple-init-params-six-raster-list.tsv"
    # method <- "analytic"
    # prev.ht <- NULL
    # maxtime <- NULL
    # outfile <- NULL
}

if (! method %in% c("analytic","CG")) { stop(usage) }

source("resistance-fns.R")
require(raster)

require(parallel)
numcores<-as.numeric(scan(pipe("cat /proc/cpuinfo | grep processor | tail -n 1 | awk '{print $3}'")))+1

if (is.null(outfile)) { outfile <- paste( subdir, "/", basename(layer.file), "-hitting-times.tsv", sep='') }

layer.names <- scan(layer.file,what="char") 

load( paste(subdir,"/",basename(layer.prefix),"_",basename(layer.file),"_","G.RData",sep='') ) # provides "G"        "update.G" "ndelta"   "ngamma"   "transfn"  "valfn"    "layers"
Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )

load( paste( subdir, "/", basename(layer.prefix), basename(layer.file), "_neighborhoods.RData", sep='' ) ) # provides 'neighborhoods'
load(paste(subdir,"/",basename(layer.prefix),"tortlocs.RData",sep='')) # provides 'locs'

# REMOVE MISSING INDIV
na.indiv <- which( is.na( locs ) )
locs <- locs[-na.indiv]
neighborhoods <- neighborhoods[-na.indiv]


##
# initial parameters?
#   in time t, 1D RW with rate r does t*r jumps, displacement has variance t*r
#   so time to move N grid sites away is sqrt(N)/r
#   so if hitting times are of order T, want r of order sqrt(N)/T

init.param.table <- read.table( param.file, header=TRUE )
init.params <- unlist( init.param.table[ match( subdir, init.param.table[,1] ), -1 ] )

G@x <- update.G(init.params)


