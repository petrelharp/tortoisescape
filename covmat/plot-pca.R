#!/usr/bin/Rscript
require(colorspace)
require(sp)
require(rgdal)
require(maps)
require(maptools)

torts <- read.csv("../1st_180_torts.csv",header=TRUE)
tort.dist.table <- read.table("../1st180_pairwise_distances_sorted_redundancy_removed.txt",header=TRUE,stringsAsFactors=FALSE)
tort.dist.table$etort1 <- factor( tort.dist.table$etort1, levels=levels(torts$EM_Tort_ID) )
tort.dist.table$etort2 <- factor( tort.dist.table$etort2, levels=levels(torts$EM_Tort_ID) )
tort.dists <- matrix( NA, nrow=nrow(torts), ncol=nrow(torts) )
tort.dists[ cbind( as.numeric( tort.dist.table$etort1 ) , as.numeric( tort.dist.table$etort2 ) ) ] <- tort.dist.table$DISTANCE
tort.dists[ cbind( as.numeric( tort.dist.table$etort2 ) , as.numeric( tort.dist.table$etort1 ) ) ] <- tort.dist.table$DISTANCE
diag(tort.dists) <- 0

# remove 12N indiv's values
torts$Northing[torts$UTM_Zone == "12N"] <- NA
torts$Easting[torts$UTM_Zone == "12N"] <- NA

rasters <- read.csv("../raster.values.tort.locations.matrix.csv",row.names=1)
stopifnot( all( rownames(rasters) == torts$EM_Tort_ID ) )
torts <- cbind(torts,rasters)

## PCA
covmat <- scan("alleleCounts100k-covmat.txt")
dim(covmat) <- c(nrow(torts),nrow(torts))

eig.covmat <- eigen(covmat)

torts$PC1 <- eig.covmat$vectors[,1]
torts$PC2 <- eig.covmat$vectors[,2]
torts$PC3 <- eig.covmat$vectors[,3]
torts$PC4 <- eig.covmat$vectors[,4]

# get elevation raster
elev.file <- "../geolayers/TIFF/10x/crop_resampled_masked_aggregated_10x_dem_30.gri"
elev <- raster(elev.file)
# and tortoise locs
load("../tort.coords.rasterGCS.Robj")
# and county lines
load("../county_tortoise_plotting_info.Robj")  # @gbradburd: how was this produced?

pdf(file="PCs-and-position.pdf", width=10, height=8, pointsize=10)
pairs( subset(torts,UTM_Zone!="12N")[c("PC1","PC2","PC3","PC4","Northing","Easting")] )
dev.off()

pdf(file="PC-maps.pdf", width=10, height=8, pointsize=10)

    cols <- diverge_hcl(64)
    pairs( torts[c("PC1","PC2","PC3","PC4")], 
        lower.panel=function(x,y,...){points(x,y,bg=cols[cut(torts$Northing,breaks=64)], pch=21, cex=2, col=grey(.20) )}, 
        upper.panel=function(x,y,...){points(x,y,bg=cols[cut(torts$Easting,breaks=64)], pch=21, cex=2, col=grey(.20) )},
        main="Northing (above), Easting (below)" )

    layout(matrix(1:4,nrow=2))
    plot( Northing ~ Easting, data=torts, bg=cols[cut(PC1,breaks=64)], main="PC1", pch=21, cex=2 )
    plot( Northing ~ Easting, data=torts, bg=cols[cut(PC2,breaks=64)], main="PC2", pch=21, cex=2 )
    plot( Northing ~ Easting, data=torts, bg=cols[cut(PC3,breaks=64)], main="PC3", pch=21, cex=2 )
    plot( Northing ~ Easting, data=torts, bg=cols[cut(PC4,breaks=64)], main="PC4", pch=21, cex=2 )

dev.off()

# check stuff for technical artifacts
for (x in grep("PC",names(torts),value=TRUE,invert=TRUE)) {
    layout(matrix(1:4,nrow=2))
    for (pc in c("PC1","PC2","PC3","PC4")) {
        plot( torts[[x]][torts$UTM_Zone!="12N"], torts[[pc]][torts$UTM_Zone!="12N"], main=paste(x,pc) )
    }
    if (is.null(locator(1))) break
}


# Does PC1 correlate with coverage?
coverages <- scan("alleleCounts500kLoci.colmeans.txt")
dim(coverages) <- c(2,nrow(torts))
torts$coverage <- colSums(coverages)
torts$mean.major <- coverages[1,]

plot( PC1 ~ coverage, data=torts )
plot( PC1 ~ mean.major, data=torts )
plot( PC1 ~ I(mean.major/coverage), data=torts )


if (FALSE) {

    ###
    # remove cluster in Ivanpah

    ivanpah <- ( torts$Location_ID %in% c("Ivanpah","ISEGS","Silver State") )

    sub.eig.covmat <- eigen( covmat[!ivanpah,!ivanpah] )
    rownames( sub.eig.covmat$vectors ) <- torts$EM_Tort_ID[ which(!ivanpah) ]

    for (k in 1:4) {
        torts[[paste("subPC",k,sep='')]] <- NA
        torts[[paste("subPC",k,sep='')]][which(!ivanpah)] <- sub.eig.covmat$vectors[,k]
    }

    pairs( subset(torts,UTM_Zone!="12N")[c("subPC1","subPC2","subPC3","subPC4","Northing","Easting")] )
}
