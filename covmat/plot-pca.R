#!/usr/bin/Rscript
require(colorspace)

torts <- read.csv("../1st_180_torts.csv",header=TRUE)
tort.dist.table <- read.table("../1st180_pairwise_distances_sorted_redundancy_removed.txt",header=TRUE,stringsAsFactors=FALSE)
tort.dist.table$etort1 <- factor( tort.dist.table$etort1, levels=levels(torts$EM_Tort_ID) )
tort.dist.table$etort2 <- factor( tort.dist.table$etort2, levels=levels(torts$EM_Tort_ID) )
tort.dists <- matrix( NA, nrow=nrow(torts), ncol=nrow(torts) )
tort.dists[ cbind( as.numeric( tort.dist.table$etort1 ) , as.numeric( tort.dist.table$etort2 ) ) ] <- tort.dist.table$DISTANCE
tort.dists[ cbind( as.numeric( tort.dist.table$etort2 ) , as.numeric( tort.dist.table$etort1 ) ) ] <- tort.dist.table$DISTANCE
diag(tort.dists) <- 0

covmat <- scan("alleleCounts100k-covmat.txt")
dim(covmat) <- c(nrow(torts),nrow(torts))

eig.covmat <- eigen(covmat)

torts$PC1 <- eig.covmat$vectors[,1]
torts$PC2 <- eig.covmat$vectors[,2]
torts$PC3 <- eig.covmat$vectors[,3]
torts$PC4 <- eig.covmat$vectors[,4]

pdf(file="PCs-and-position.pdf", width=10, height=8, pointsize=10)
pairs( subset(torts,UTM_Zone!="12N")[c("PC1","PC2","PC3","PC4","Northing","Easting")] )
dev.off()

if (FALSE) {

    layout(t(1:2))
    plot( PC2 ~ PC1, data=torts, col=terrain_hcl(64)[cut(Northing,breaks=64)], main="Northing", subset=UTM_Zone!="12N", pch=20, cex=2 )
    plot( PC2 ~ PC1, data=torts, col=terrain_hcl(64)[cut(Easting,breaks=64)], main="Easting", subset=UTM_Zone!="12N", pch=20, cex=2 )

    layout(t(1:2))
    plot( Northing ~ Easting, data=torts, col=terrain_hcl(64)[cut(PC1,breaks=64)], main="PC1", subset=UTM_Zone!="12N", pch=20, cex=2 )
    plot( Northing ~ Easting, data=torts, col=terrain_hcl(64)[cut(PC2,breaks=64)], main="PC2", subset=UTM_Zone!="12N", pch=20, cex=2 )

    plot( Northing ~ Easting, data=torts, subset=UTM_Zone!="12N", pch=20, cex=2, col=Location_ID )
}

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
