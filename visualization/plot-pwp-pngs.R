#!/usr/bin/Rscript
library(methods)  # required for get( ) below... ???

usage <- "Makes many images for a single data frame of pairwise comparisons: one for each thing being compared.
Usage:
    Rscript plot-pwp-pngs.R (distance file 1) (distance file 2) (sample info directory) (output directory)
"

argvec <- if (interactive()) { scan(what="char") } else { commandArgs(TRUE) }

if (length(argvec)<4) { stop(usage) }

dist1.file <- argvec[1]
dist2.file <- argvec[2]
indir <- argvec[3]
outdir <- argvec[4]

dist1 <- read.csv(dist1.file,header=TRUE,stringsAsFactors=FALSE)
dist2 <- read.csv(dist2.file,header=TRUE,stringsAsFactors=FALSE)
dir.create(outdir,showWarnings=FALSE)

# locations
coord.obj <- load(file.path(indir,"geog_coords.RData"))
coords <- get(coord.obj)
tort.ids <- row.names(coords)

# read in other info
pcs <- read.csv(file.path(indir,"pcs.csv"),header=TRUE,stringsAsFactors=FALSE)
stopifnot( all( tort.ids %in% pcs$etort ) )
pc.cols <- adjustcolor( ifelse( pcs$PC1[match(tort.ids,pcs$etort)] > 0, "blue", "purple" ), .75 )

# reorder so everything matches

require(raster)
layer <- raster("../visualization/crm_dem_30")

for (tid in tort.ids) {
    cat(tid,"\n")
  png( file=file.path(outdir,paste(tid,".png",sep='')), width=8*144, height=4*144, pointsize=10, res=144 )
    usethese <- ( dist1$etort1 == tid ) | ( dist1$etort2 == tid )
    otherone <- ifelse( dist1$etort1[usethese] == tid, dist1$etort2[usethese], dist1$etort1[usethese] )
    thiscolors <- pc.cols[ match(otherone,pcs$etort) ]
    layout(t(1:2))
    plot(layer)
    points(coords,pch=20,cex=1,col=pc.cols)
    points(coords[match(tid,tort.ids)],cex=2,col='red')
    plot( dist1[,3], dist2[,3], pch=20, cex=.5, 
        col=adjustcolor("black",.25), xlab=dist1.file, ylab=dist2.file )
    points( dist1[,3][usethese], dist2[,3][usethese], pch=20, col=thiscolors, cex=1.5 )
  dev.off()
}


