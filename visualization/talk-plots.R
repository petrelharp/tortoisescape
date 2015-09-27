#!/usr/bin/Rscript

dist1.file <- "../tort_272_info/geog_distance.csv"
dist2.file <- "../tort_272_info/all_angsd_snps.pwp.csv"
indir <- "../tort_272_info"

dist1 <- read.csv(dist1.file,header=TRUE,stringsAsFactors=FALSE)
dist2 <- read.csv(dist2.file,header=TRUE,stringsAsFactors=FALSE)

# locations
coord.obj <- load(file.path(indir,"geog_coords.RData"))
coords <- get(coord.obj)
tort.ids <- row.names(coords)

# read in other info
pcs <- read.csv(file.path(indir,"pcs.csv"),header=TRUE,stringsAsFactors=FALSE)
stopifnot( all( tort.ids %in% pcs$etort ) )
base.cols <- adjustcolor( c("blue","purple","red"), 0.25 )  # north, south, between
pc.cols <- base.cols[ ifelse( pcs$PC1[match(tort.ids,pcs$etort)] > 0, 1, 2 ) ]

relatives <- ( (dist2$etort1==dist2$etort2) | ( dist2[,3] < quantile(subset(dist2,etort1==etort2)[,3],0.75) ) )

require(raster)
layer <- raster("../visualization/dem_30")

player <- function (main='') { plot(layer,legend=FALSE,xlab="",ylab="",xaxt="n",yaxt="n",legend.mar=0,box=FALSE,main=main) }

pdf(file="everyone-pwp.pdf",width=2.5,height=5,pointsize=10)
    usethese <- !relatives
    north1 <- ( pcs$PC1[match(dist2$etort1[usethese],pcs$etort)] > 0 )
    north2 <- ( pcs$PC1[match(dist2$etort2[usethese],pcs$etort)] > 0 )
    distcolors <- base.cols[ ifelse( north1&north2, 1, ifelse( (!north1)&(!north2), 2, 3 ) ) ]
    layout((1:2))
    par(mar=c(0,0,0.5,0))
    player()
    points(coords,pch=20,col=pc.cols)
    par(mar=c(2.5,2.5,0.5,0.5))
    plot( dist1[usethese,3], dist2[usethese,3], pch=20, cex=.25, 
       col=distcolors, xlab="geog dist (km)", ylab="divergence",
       mgp=c(2.0,0.75,0) )
dev.off()

for (tid in c(1, 100, 200) {
  png( file=file.path(outdir,paste(tour.labels[k], "_",gsub("[^0-9a-z-]","_",tid),".png",sep='')), width=12*144, height=4*144, pointsize=10, res=144 )
    layout(t(1:3))
    opar <- par(mar=c(1,1,2,1))
    usethese <- ( dist2$etort1 != dist2$etort2 ) & ( ( dist2$etort1 == tid ) | ( dist2$etort2 == tid ) )
    otherone <- ifelse( dist2$etort1[usethese] == tid, dist2$etort2[usethese], dist2$etort1[usethese] )
    thiscolors <- pc.cols[ match(otherone,tort.ids) ]
    player(paste(tid," similarities"))
    points(coords[match(otherone,tort.ids)],pch=20,cex=sfn(dist2[,3][usethese]),col=thiscolors)
    points(coords[match(tid,tort.ids)],pch="*",cex=2)
    player(paste(tid," distances"))
    points(coords[match(otherone,tort.ids)],pch=20,cex=dfn(dist2[,3][usethese]),col=thiscolors)
    points(coords[match(tid,tort.ids)],pch="*",cex=2)
    par(opar)
    plot( dist1[,3], dist2[,3], pch=20, cex=.5, 
        col=adjustcolor("black",.25), xlab=dist1.file, ylab=dist2.file )
    points( dist1[,3][usethese], dist2[,3][usethese], pch=20, col=thiscolors, cex=1.5 )
  dev.off()
}
