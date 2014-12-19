#!/usr/bin/Rscript

usage <- "
Compute hitting times at a given set of parameters, add noise, and write out the pairwise times for sample locations (to mimic the data).
Usage:
    Rscript sim-hitting-times.R (layer prefix) (subdir) (layer file) (parameter file) (size of noise, as a proportion) (output file)
e.g.
    Rscript sim-hitting-times.R ../geolayers/multigrid/256x/crm_ 256x  six-raster-list test01/six-params.tsv 0.05 test01/256x/six-raster-list-sim-hts.tsv

"

if (!interactive()) {
    if (length(commandArgs(TRUE))<6) { stop(usage) }
    layer.prefix <- commandArgs(TRUE)[1]
    subdir <- commandArgs(TRUE)[2]
    layer.file <- commandArgs(TRUE)[3]
    param.file <- commandArgs(TRUE)[4]
    noise.fac <- as.numeric(commandArgs(TRUE)[5])
    output.file <- commandArgs(TRUE)[6]
} else {
    layer.prefix <- c("../geolayers/multigrid/256x/crm_")
    subdir <- "256x"
    layer.file <- "six-raster-list"
    param.file <- "test01/six-params.tsv"
    noise.fac <- 0.05
    output.file <- "test01/256x/six-raster-list-sim-hts.tsv"
}
cat("Parameters: \n")
invisible( lapply( c("layer.prefix","subdir","layer.file","param.file","noise.fac","output.file"), function (x) { cat("  ", x, " : ", get(x), "\n") } ) )
cat("\n")
layer.names <- scan(layer.file,what="char") 

load(paste(subdir,"/",basename(layer.prefix),basename(layer.file),"-","setup.RData",sep=''))

source("resistance-fns.R")
require(raster)

require(parallel)
numcores<-getcores()

init.params.table <- read.table(param.file,header=TRUE)
init.params <- as.numeric( init.params.table[match(subdir,init.params.table[,1]),-1] )
names(init.params) <- colnames(init.params.table)[-1]

G@x <- update.G(init.params)

hts <- hitting.analytic( neighborhoods, G-diag(rowSums(G)), numcores=getcores() )

## DO NOT SYMMETRIZE
# sym.hts <- ( hts[locs,] + t(hts[locs,]) ) / 2
sym.hts <- hts[locs,]
noisy.df <- sym.hts.df <- data.frame(
        etort1=row(sym.hts)[upper.tri(sym.hts,diag=TRUE)],
        etort2=col(sym.hts)[upper.tri(sym.hts,diag=TRUE)],
        DISTANCE=sym.hts[upper.tri(sym.hts,diag=TRUE)]
    )

noisy.df$DISTANCE <- sym.hts.df$DISTANCE * exp( rnorm(nrow(sym.hts.df),sd=noise.fac) )

dir.create( dirname(output.file), recursive=TRUE, showWarnings=FALSE )
cat("Writing out to ", output.file, " .\n")
write.table( noisy.df, file=output.file, row.names=FALSE )
