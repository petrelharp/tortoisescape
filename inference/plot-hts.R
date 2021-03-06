#!/usr/bin/Rscript

usage <- '
    Read in and plot fitted hitting times.  Usage:
        Rscript plot-hts.R (layer prefix) (subdirectory) (layer file name) (tsv of hitting times) (base name for output pngs)
    e.g.
        Rscript plot-hts.R ../geolayers/multigrid/64x/crm_ 64x dem-layer-list 64x/dem-layer-list-hitting-times.tsv dem/64x_hts_

'

if (!interactive()) {
    if (length(commandArgs(TRUE))<5) { cat(usage); q() }
    layer.prefix <- commandArgs(TRUE)[1]
    subdir <- commandArgs(TRUE)[2]
    layer.file <- commandArgs(TRUE)[3]
    ht.file <- commandArgs(TRUE)[4]
    outbase <-  commandArgs(TRUE)[5]
} else {
    layer.prefix <- "../geolayers/multigrid/64x/crm_"
    subdir <- "64x"
    layer.file <- "dem-layer-list"
    ht.file <- "64x/dem-layer-list-hitting-times.tsv"
    outbase <- "dem/64x_hts_"
}

source("resistance-fns.R")
require(raster)

require(parallel)
numcores <- getcores()

load( paste( subdir, "/", basename(layer.prefix), basename(layer.file), "_neighborhoods.RData", sep='' ) ) # provides 'neighborhoods'

# layer.names <- scan(layer.file,what="char") 
# load( paste(subdir,"/",basename(layer.prefix),"_",basename(layer.file),"_","G.RData",sep='') ) # provides "G"    "Gjj"    "update.G" "ndelta"   "ngamma"   "transfn"  "valfn"    "layers"

# HAVE REMOVED MISSING INDIV
load(paste(subdir,"/",basename(layer.prefix),basename(layer.file),"_tortlocs.RData",sep='')) # provides 'locs'
orig.locs <- locs
na.indiv <- which( is.na( orig.locs ) )
locs <- orig.locs[-na.indiv]
neighborhoods <- lapply(neighborhoods[-na.indiv],function (x) { x[!is.na(x)] })
tort.nums <- seq_along(orig.locs)[-na.indiv]

load( paste(subdir, "/", basename(layer.prefix),"_", basename(layer.file),"_nonmissing.RData",sep='') ) # provides nonmissing
ph <- plot.ht.fn(layer.prefix,nonmissing,"annual_precip",par.args=list(mar=c(2,2,3,2)+.1))

hts <- as.matrix( read.table(ht.file,header=TRUE) )

dir.create(file.path(subdir,dirname(outbase)),recursive=TRUE)

mclapply( 1:ncol(hts), function (k) {
        # 'locs' are indices of tortoise locations
        ymax <- 1.5*max(hts[locs,k])  # 1.5 times maximum hitting time to another tortoise location
        yy <- hts[,k]
        yy[yy>ymax] <- NA
        outfile <- file.path(subdir,paste(outbase,"_",tort.nums[k],".png",sep=""))
        png( file=outfile, res=144, width=5*144, height=4*144, pointsize=10 )
        ph( yy, main=paste("Tortoise",tort.nums[k]) )
        dev.off()
        return(outfile)
    }, mc.cores=numcores )
