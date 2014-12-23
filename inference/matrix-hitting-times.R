#!/usr/bin/Rscript

usage <- '
    Read in fitted hitting times and extract the matrix of hitting times by sample. Usage:
        Rscript matrix-hitting-times.R (layer prefix) (subdirectory) (layer file name) (tsv of hitting times) (output file name)
    e.g.
        Rscript matrix-hitting-times.R ../geolayers/multigrid/128x/crm_ 128x six-raster-list 128x/six-raster-list-hitting-times.tsv 128x/six-raster-list-hitting-times-torts.tsv

'

if (!interactive()) {
    if (length(commandArgs(TRUE))<5) { cat(usage); q() }
    layer.prefix <- commandArgs(TRUE)[1]
    subdir <- commandArgs(TRUE)[2]
    layer.file <- commandArgs(TRUE)[3]
    ht.file <- commandArgs(TRUE)[4]
    outfile <-  commandArgs(TRUE)[5]
} else {
    layer.prefix <- "../geolayers/multigrid/128x/crm_"
    subdir <- "128x"
    layer.file <- "six-raster-list"
    ht.file <- "128x/six-raster-list-hitting-times.tsv"
    outfile <- "128x/six-raster-list-hitting-times-torts.tsv"
}

load(paste(subdir,"/",basename(layer.prefix),basename(layer.file),"_tortlocs.RData",sep='')) # provides 'locs'
orig.locs <- locs
na.indiv <- which( is.na( orig.locs ) )
locs <- orig.locs[-na.indiv]
tort.nums <- seq_along(orig.locs)[-na.indiv]

require(raster)
load("../tort.coords.rasterGCS.Robj")
all.tort.ids <- rownames(tort.coords.rasterGCS@coords)
use.tort.ids <- all.tort.ids[-na.indiv]

hts <- as.matrix( read.table(ht.file,header=TRUE) )

out.hts <- as.data.frame( t(combn(all.tort.ids,2)), stringsAsFactors=FALSE )
names(out.hts) <- c("etort1","etort2")
out.hts$ht <- hts[ cbind( locs[match(out.hts$etort1,use.tort.ids)], match(out.hts$etort2,use.tort.ids) ) ]

write.table(out.hts,file=outfile,row.names=FALSE)
