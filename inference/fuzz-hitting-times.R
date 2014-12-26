#!/usr/bin/Rscript

usage <- "
Read in hitting times at a given set of parameters, add noise, and write out the pairwise times for sample locations (to mimic the data).
Usage:
    Rscript sim-hitting-times.R (layer prefix) (subdir) (layer file) (input file) (size of noise, as a proportion) (output file)
e.g.
    Rscript sim-hitting-times.R ../geolayers/multigrid/256x/crm_ 256x  six-raster-list test_six_layers/256x/six-raster-list-hitting-times-full.tsv 0.01 test_six_layers/256x/six-raster-list-sim-0_01-hts.tsv

"

argvec <- if (!interactive()) { commandArgs(TRUE) } else { scan(what='char') }
if (length(argvec)<6) { stop(usage) }
layer.prefix <- argvec[1]
subdir <- argvec[2]
layer.file <- argvec[3]
input.file <- argvec[4]
noise.fac <- as.numeric(argvec[5])
output.file <- argvec[6]

cat("Parameters: \n")
invisible( lapply( c("layer.prefix","subdir","layer.file","input.file","noise.fac","output.file"), function (x) { cat("  ", x, " : ", get(x), "\n") } ) )
cat("\n")

layer.names <- scan(layer.file,what="char") 

load(paste(subdir,"/",basename(layer.prefix),basename(layer.file),"-","setup.RData",sep=''))

source("resistance-fns.R")
require(raster)

require(parallel)
numcores<-getcores()

# read in hitting times as output by make-resistance-distances
hts <- read.full.hts( input.file, locs )

## DO NOT SYMMETRIZE
# sym.hts <- ( hts[locs,] + t(hts[locs,]) ) / 2
sym.hts <- hts[locs,]
noisy.hts <- sym.hts * exp( rnorm(length(sym.hts))*noise.fac )

dir.create( dirname(output.file), recursive=TRUE, showWarnings=FALSE )
write.sub.hts( noisy.hts, file=output.file )
