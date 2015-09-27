#!/usr/bin/Rscript

dist1.file <- "../tort_272_info/geog_distance.csv"
dist2.file <- "../tort_272_info/all_angsd_snps.pwp.csv"
indir <- "../tort_272_info"

dist1 <- read.csv(dist1.file,header=TRUE,stringsAsFactors=FALSE)
dist2 <- read.csv(dist2.file,header=TRUE,stringsAsFactors=FALSE)

load("../visualization/county_lines.Robj")  # provides county_lines

# locations
coord.obj <- load(file.path(indir,"geog_coords.RData"))
coords <- get(coord.obj)
tort.ids <- row.names(coords)

# read in other info
pcs <- read.csv(file.path(indir,"pcs.csv"),header=TRUE,stringsAsFactors=FALSE)
stopifnot( all( tort.ids %in% pcs$etort ) )
base.cols <- c("blue","purple","red")  # north, south, between
pc.cols <- adjustcolor(base.cols,0.75)[ ifelse( pcs$PC1[match(tort.ids,pcs$etort)] > 0, 1, 2 ) ]

relatives <- ( (dist2$etort1==dist2$etort2) | ( dist2[,3] < quantile(subset(dist2,etort1==etort2)[,3],0.75) ) )

require(raster)
layer <- raster("../visualization/dem_30")

player <- function (main='') { 
    plot(layer,legend=FALSE,xlab="",ylab="",xaxt="n",yaxt="n",legend.mar=0,box=FALSE,main=main) 
    lines(county_lines,lwd=0.5)
}

pdf(file="everyone-pwp.pdf",width=5.5,height=5.5/1.7,pointsize=10)
    usethese <- !relatives
    north1 <- ( pcs$PC1[match(dist2$etort1[usethese],pcs$etort)] > 0 )
    north2 <- ( pcs$PC1[match(dist2$etort2[usethese],pcs$etort)] > 0 )
    layout(t(1:2))
    par(mar=c(2.5,2.5,0.5,0.5))
    player()
    points(coords,pch=20,col=pc.cols)
    plot( dist1[usethese,3], dist2[usethese,3], pch=20, cex=.25, 
       col=adjustcolor("black",0.25), xlab="geog dist (km)", ylab="divergence",
       mgp=c(1.6,0.75,0) )
dev.off()

pdf(file="everyone-pwp-vertical.pdf",width=2.5,height=5,pointsize=10)
    usethese <- !relatives
    north1 <- ( pcs$PC1[match(dist2$etort1[usethese],pcs$etort)] > 0 )
    north2 <- ( pcs$PC1[match(dist2$etort2[usethese],pcs$etort)] > 0 )
    distcolors <- adjustcolor(base.cols,0.25)[ ifelse( north1&north2, 1, ifelse( (!north1)&(!north2), 2, 3 ) ) ]
    layout((1:2))
    par(mar=c(0,0,0.5,0))
    player()
    points(coords,pch=20,col=pc.cols)
    par(mar=c(2.5,2.5,0.5,0.5))
    plot( dist1[usethese,3], dist2[usethese,3], pch=20, cex=.25, 
       col=distcolors, xlab="geog dist (km)", ylab="divergence",
       mgp=c(1.6,0.75,0) )
dev.off()

mindist <- min(dist2[!relatives,3])
sddist <- sd(dist2[!relatives,3])
sfn <- function (x,max.cex=7) {
    max.cex/( 1 + exp( (x-mindist)/sddist ) )
}

for (tid in paste("etort-",c(285,240,35,273,57,229,191),sep='')) {
  pdf( file=paste("pwp_",tid,".pdf",sep=''), width=5.5, height=5.5/1.7, pointsize=10 )
    layout(t(1:2))
    par(mar=c(2.5,2.5,0.5,0.5))
    usethese <- ( dist2$etort1 != dist2$etort2 ) & ( ( dist2$etort1 == tid ) | ( dist2$etort2 == tid ) )
    otherone <- ifelse( dist2$etort1[usethese] == tid, dist2$etort2[usethese], dist2$etort1[usethese] )
    thiscolors <- pc.cols[ match(otherone,tort.ids) ]
    player()
    points(coords[match(otherone,tort.ids)],pch=20,cex=sfn(dist2[,3][usethese]),col=thiscolors)
    points(coords[match(tid,tort.ids)],pch="*",cex=2)
    plot( dist1[!relatives,3], dist2[!relatives,3], pch=20, cex=.5, 
        col=adjustcolor("black",.25),
       xlab="geog dist (km)", ylab="divergence",
       mgp=c(1.6,0.75,0) )
    points( dist1[,3][usethese], dist2[,3][usethese], pch=20, col=thiscolors, cex=1.5 )
  dev.off()
}
