#!/usr/bin/Rscript
library(methods)  # required for get( ) below... ???

usage <- "Makes many images for a single data frame of pairwise comparisons: one for each thing being compared.
Usage:
    Rscript distance-maps.R (distance file) (sample info directory) (output directory)
"

argvec <- if (interactive()) { scan(what="char") } else { commandArgs(TRUE) }

if (length(argvec)<3) { stop(usage) }

dist.file <- argvec[1]
indir <- argvec[2]
outdir <- argvec[3]

dist <- read.csv(dist.file,header=TRUE,stringsAsFactors=FALSE)
dir.create(outdir,showWarnings=FALSE)

# locations
coord.obj <- load(file.path(indir,"geog_coords.RData"))
coords <- get(coord.obj)
tort.ids <- row.names(coords)

# read in other info
pcs <- read.csv(file.path(indir,"pcs.csv"),header=TRUE,stringsAsFactors=FALSE)
stopifnot( all( tort.ids %in% pcs$etort ) )
pc.cols <- adjustcolor( ifelse( pcs$PC1[match(tort.ids,pcs$etort)] > 0, "blue", "purple" ), .75 )


require(raster)
layer <- raster("../visualization/dem_30")

player <- function (main='') { plot(layer,legend=FALSE,xlab="",ylab="",xaxt="n",yaxt="n",legend.mar=0,box=FALSE,main=main) }

##
# find a good ordering

require(TSP)
xy <- coordinates(coords)
etsp <- ETSP( xy, labels=rownames(xy) )
tour <- solve_TSP( etsp, method="linkern" )
tour.labels <- t(outer(letters,letters,paste,sep=''))[seq_len(length(tour))]

png( file=file.path(outdir,"plot-order.png"), width=4*144, height=4*144, pointsize=10, res=144 )
player("tour")
segments(x0=xy[tour,1],x1=xy[c(tour[-1],tour[1]),1],
    y0=xy[tour,2],y1=xy[c(tour[-1],tour[1]),2])
dev.off()



# First plot self-comparisons
png( file=file.path(outdir,"self-comparisons.png"), width=8*144, height=4*144, pointsize=10, res=144 )
layout(t(1:2))
par(mar=c(1,1,2,1))
usethese <- ( dist$etort1 == dist$etort2 )
thiscolors <- pc.cols[ match(dist$etort1,tort.ids) ]
x <- dist[usethese,3]
player("self-similarities")
points(coords,pch=20,col=pc.cols,cex=3/(1+exp((x-min(x))/sd(x))))
player("self-distances")
points(coords,pch=20,col=pc.cols,cex=(x-min(x))/(2*sd(x)))
dev.off()

relatives <- ( (dist$etort1==dist$etort2) | ( dist[,3] < quantile(subset(dist,etort1==etort2)[,3],0.75) ) )
mindist <- min(dist[!relatives,3])
sddist <- sd(dist[!relatives,3])
sfn <- function (x,max.cex=7) {
    max.cex/( 1 + exp( (x-mindist)/sddist ) )
}
dfn <- function (x) {
    (x-mindist)/(2*sddist)
}


for (k in seq_along(tort.ids)) {
    tid <- tort.ids[tour[k]]
    cat(tid,"\n")
  png( file=file.path(outdir,paste(tour.labels[k], "_",gsub("[^0-9a-z-]","_",tid),".png",sep='')), width=8*144, height=4*144, pointsize=10, res=144 )
    layout(t(1:2))
    par(mar=c(1,1,2,1))
    usethese <- ( dist$etort1 != dist$etort2 ) & ( ( dist$etort1 == tid ) | ( dist$etort2 == tid ) )
    otherone <- ifelse( dist$etort1[usethese] == tid, dist$etort2[usethese], dist$etort1[usethese] )
    thiscolors <- pc.cols[ match(otherone,tort.ids) ]
    player(paste(tid," similarities"))
    points(coords[match(otherone,pcs$etort)],pch=20,cex=sfn(dist[,3][usethese]),col=thiscolors)
    points(coords[match(tid,tort.ids)],pch="*",cex=2)
    player(paste(tid," distances"))
    points(coords[match(otherone,tort.ids)],pch=20,cex=dfn(dist[,3][usethese]),col=thiscolors)
    points(coords[match(tid,tort.ids)],pch="*",cex=2)
  dev.off()
}
