#!/usr/bin/Rscript
library(methods)  # required for get( ) below... ???

usage <- "Makes many images for a single data frame of pairwise comparisons: one for each thing being compared.
Usage:
    Rscript ms-maps.R (directory with simulation output in)
"

argvec <- if (interactive()) { scan(what="char") } else { commandArgs(TRUE) }

if (length(argvec)<1) { stop(usage) }

indir <- argvec[1]
dist1.file <- file.path(indir,"mean-divergence.csv")
dist2.file <- file.path(indir,"geog-distance.csv")
outdir <- file.path(indir,"pngs")

dist1 <- read.csv(dist1.file,header=TRUE,stringsAsFactors=FALSE)
dist2 <- read.csv(dist2.file,header=TRUE,stringsAsFactors=FALSE)
dir.create(outdir,showWarnings=FALSE)

# locations
coord.obj <- load(file.path(indir,"geog_coords.RData"))
coords <- get(coord.obj)
tort.ids <- row.names(coords)

# read in other info
pcs <- read.csv(file.path(indir,"pcs.csv"),header=TRUE,stringsAsFactors=FALSE)
stopifnot( all( tort.ids %in% pcs[,1] ) )
pc.cols <- adjustcolor( ifelse( pcs[,2][match(tort.ids,pcs[,1])] > 0, "blue", "purple" ), .75 )

require(raster)
layer <- crop( raster("../visualization/nussear_masked.grd"), coords )

player <- function (main='') { plot(layer,legend=FALSE,xlab="",ylab="",xaxt="n",yaxt="n",legend.mar=0,box=FALSE,main=main) }

##
# find a good ordering

require(TSP)
xy <- jitter(coordinates(coords))  # doesn't like it if points coincide
etsp <- ETSP( xy, labels=rownames(xy) )
tour <- solve_TSP( etsp, method="linkern" )
tour.labels <- t(outer(letters,letters,paste,sep=''))[seq_len(length(tour))]


mindist <- min(dist2[,3])
sddist <- sd(dist2[,3])
sfn <- function (x,max.cex=7) {
    max.cex/( 1 + exp( (x-mindist)/sddist ) )
}


for (k in seq_along(tort.ids)) {
    tid <- tort.ids[tour[k]]
    cat(tid,"\n")
  png( file=file.path(outdir,paste(tour.labels[k], "_",gsub("[^0-9a-z-]","_",tid),".png",sep='')), width=12*144, height=4*144, pointsize=10, res=144 )
    layout(t(1:2))
    opar <- par(mar=c(1,1,2,1))
    usethese <- ( dist2[,1] != dist2[,2] ) & ( ( dist2[,1] == tid ) | ( dist2[,2] == tid ) )
    otherone <- ifelse( dist2[,1][usethese] == tid, dist2[,2][usethese], dist2[,1][usethese] )
    thiscolors <- pc.cols[ match(otherone,tort.ids) ]
    player(paste(tid," similarities"))
    points(coords[match(otherone,tort.ids)],pch=20,cex=sfn(dist2[,3][usethese]),col=thiscolors)
    points(coords[match(tid,tort.ids)],pch="*",cex=2)
    par(opar)
    plot( dist1[,3], dist2[,3], pch=20, cex=.5, 
        col=adjustcolor("black",.25), xlab=dist1.file, ylab=dist2.file )
    points( dist1[,3][usethese], dist2[,3][usethese], pch=20, col=thiscolors, cex=1.5 )
  dev.off()
}
