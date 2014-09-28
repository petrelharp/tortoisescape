#!/usr/bin/Rscript

load("tort.coords.rasterGCS.Robj")

require(raster)
layer <- raster("geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_dem_30.gri")

torts <- read.csv("1st_180_torts.csv",header=TRUE,stringsAsFactors=FALSE)
torts$EM_Tort_ID <- factor( torts$EM_Tort_ID , levels=torts$EM_Tort_ID )
nind <- nrow(torts)

dists <- read.csv("pairwise-normalized-pi.csv",header=TRUE,stringsAsFactors=FALSE) # had DISTANCE, pi, and npi
dists$etort1 <- factor( dists$etort1 , levels=torts$EM_Tort_ID )
dists$etort2 <- factor( dists$etort2 , levels=torts$EM_Tort_ID )

pcs <- read.csv("covmat/tort-PCs.csv",stringsAsFactors=FALSE)
pcs$X <- factor( pcs$X, levels=torts$EM_Tort_ID )
pc.cols <- adjustcolor( ifelse( pcs$PC1 > 0, "blue", "purple" ), .75 )
stopifnot( all(pcs$X==torts$EM_Tort_ID) )

for (k in 1:nind) { 
    png( file=paste("pngs/",torts$EM_Tort_ID[k],"-ibd.png",sep=''), width=12*144, height=4*144, pointsize=10, res=144 )
    usethese <- ( dists$etort1 == torts$EM_Tort_ID[k] ) | ( dists$etort2 == torts$EM_Tort_ID[k] )
    otherone <- torts$EM_Tort_ID[ifelse( dists$etort1[usethese] == torts$EM_Tort_ID[k], dists$etort2[usethese], dists$etort1[usethese] )]
    thiscolors <- pc.cols[ match(otherone,pcs$X) ]
    layout(t(1:3))
    plot(layer)
    points(tort.coords.rasterGCS,pch=20,cex=1,col=pc.cols)
    points(tort.coords.rasterGCS[k],cex=2,col='red')
    plot( dists$DISTANCE, dists$pi, pch=20, cex=.5, 
        col=adjustcolor("black",.25), xlab="geographic distance (km)", ylab="raw pairwise divergence" )
    points( dists$DISTANCE[usethese], dists$pi[usethese], pch=20, col=thiscolors, cex=1.5 )
    plot( dists$DISTANCE, dists$npi, pch=20, cex=.5, 
        col=adjustcolor("black",.25), xlab="geographic distance (km)", ylab="adjusted pairwise divergence" )
    points( dists$DISTANCE[usethese], dists$npi[usethese], pch=20, col=thiscolors, cex=1.5 )
    # if (is.null(locator(1))) { break }
    dev.off()
}
