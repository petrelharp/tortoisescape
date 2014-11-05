#!/usr/bin/Rscript

##
# convert hitting times on one grid to another


source("resistance-fns.R")
require(raster)
rasterOptions(tmpdir=".")


require(parallel)
numcores<-getcores()


if (!interactive()) {
    layer.prefix.1 <- commandArgs(TRUE)[1]
    layer.prefix.2 <- commandArgs(TRUE)[2]
    subdir.1 <- commandArgs(TRUE)[3]
    subdir.2 <- commandArgs(TRUE)[4]
    layer.file <- commandArgs(TRUE)[5]
    prev.ht <- commandArgs(TRUE)[6]
    ag.fact <- as.numeric( commandArgs(TRUE)[7] )
} else {
    layer.prefix.1 <- "../geolayers/TIFF/500x/500x_"
    layer.prefix.2 <- "../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_"
    subdir.1 <- "500x"
    subdir.2 <- "100x"
    layer.file <- "six-raster-list"
    prev.ht <- "500x/six-raster-list-hitting-times.tsv"
    ag.fact <- 5
}

hts <- read.table(prev.ht,header=TRUE)

env <- new.env()
load( paste(subdir.1, "/", basename(layer.prefix.1), "_", basename(layer.file), "_nonmissing.RData",sep=''), envir=env ) # provides nonmissing
nonmissing.1 <- with( env, nonmissing )
load( paste(subdir.2, "/", basename(layer.prefix.2), "_", basename(layer.file), "_nonmissing.RData",sep=''), envir=env ) # provides nonmissing
nonmissing.2 <- with( env, nonmissing )

layer.name <- "dem_30"
layer.1 <- raster(paste(layer.prefix.1,layer.name,sep=''))
layer.2 <- raster(paste(layer.prefix.2,layer.name,sep=''))

values(layer.1)[-nonmissing.1] <- NA
values(layer.2)[-nonmissing.2] <- NA

##

# omit aggregation error check for larger grids
checkit <- ( subdir.1 == "500x" )

# do the disaggregation
new.hts <- upsample.hts( hts, ag.fact, layer.1, nonmissing.1, layer.2, nonmissing.2, checkit, numcores=numcores )

outfile <- paste( subdir.2, "/", basename(subdir.1), "-aggregated-hitting-times.tsv", sep='')
write.table( new.hts, file=outfile, row.names=FALSE )
cat("Writing output to ", outfile, " .\n")

if (FALSE) {
    layout(t(1:4))
    plot(layer.1)
    plot(layer.1.dis)
    plot(layer.1.dis.res)
    plot(layer.1.dis.res-layer.1.dis)
    # plot(layer.2)
}
