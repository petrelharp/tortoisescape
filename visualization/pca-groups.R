#!/usr/bin/Rscript

load("../tort_272_info/tort.coords.rasterGCS.Robj")

require(raster)
layer <- raster("../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_dem_30.gri")

torts <- read.csv("../tort_272_info/1st_180_torts.csv",header=TRUE,stringsAsFactors=FALSE)
torts$EM_Tort_ID <- factor( torts$EM_Tort_ID , levels=torts$EM_Tort_ID )
nind <- nrow(torts)

dists <- read.csv("../pairwise-normalized-pi.csv",header=TRUE,stringsAsFactors=FALSE) # had DISTANCE, pi, and npi
dists$etort1 <- factor( dists$etort1 , levels=torts$EM_Tort_ID )
dists$etort2 <- factor( dists$etort2 , levels=torts$EM_Tort_ID )

pcs <- read.csv("../covmat/tort-PCs.csv",stringsAsFactors=FALSE)
names(pcs)[1] <- "etort"
pcs$etort <- factor( pcs$etort, levels=torts$EM_Tort_ID )
stopifnot( all(pcs$etort==torts$EM_Tort_ID) )

if (interactive()) pairs(pcs[,-1])

kmean.list <- lapply( 2:9, function (nclusts) {
    kmeans( x=pcs[,paste("PC",1:5,sep='')], centers=nclusts )
    } )

if (FALSE) {
layout( matrix( seq_along(kmean.list), nrow=2 ) )
for (k  in seq_along(kmean.list)) {
    plot(layer)
    points( tort.coords.rasterGCS, pch=20, cex=1.5, col=adjustcolor(kmean.list[[k]]$cluster,.75) )
}
}

pdf(file="plots/pca-groups.pdf",width=12,height=4,pointsize=10)
layout(matrix(1:6,nrow=2))
for (k  in seq_along(kmean.list)) {
    plot(layer)
    cols <- adjustcolor(kmean.list[[k]]$cluster,.75)
    points( tort.coords.rasterGCS, pch=20, cex=1.5, col=cols )
    plot( pcs$PC1, pcs$PC2, col=cols, pch=20 )
    plot( pcs$PC3, pcs$PC2, col=cols, pch=20 )
    plot( pcs$PC3, pcs$PC1, col=cols, pch=20 )
    plot( pcs$PC4, pcs$PC2, col=cols, pch=20 )
    plot( pcs$PC5, pcs$PC2, col=cols, pch=20 )
}
dev.off()
