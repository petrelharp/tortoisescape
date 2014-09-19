#!/usr/bin/Rscript
require(colorspace)
require(sp)
require(rgdal)
require(maps)
require(maptools)
require(raster)

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
# covmat <- scan("../tortGen/exampleOutput/alleleCounts100k-covmat.txt")
covmat <- scan("../tortGen/exampleOutput/alleleCounts500kLoci-covmat.txt")
dim(covmat) <- c(nrow(torts),nrow(torts))

eig.covmat <- eigen(covmat)

# proportion of variance explained
eig.covmat$values / sum(eig.covmat$values)
cat( paste(formatC( eig.covmat$values / sum(eig.covmat$values)*100, digits=3 )[1:8],"%",sep=''), "\n" )

torts$PC1 <- eig.covmat$vectors[,1]
torts$PC2 <- eig.covmat$vectors[,2]
torts$PC3 <- eig.covmat$vectors[,3]
torts$PC4 <- eig.covmat$vectors[,4]
torts$PC5 <- eig.covmat$vectors[,5]
torts$PC6 <- eig.covmat$vectors[,6]
torts$PC7 <- eig.covmat$vectors[,7]
torts$PC8 <- eig.covmat$vectors[,8]

# get elevation raster
elev.file <- "../geolayers/TIFF/10x/crop_resampled_masked_aggregated_10x_dem_30.gri"
elev <- raster(elev.file)
# and tortoise locs
load("../tort.coords.rasterGCS.Robj")
stopifnot( all( row.names(tort.coords.rasterGCS) == torts$EM_Tort_ID ) )
# and county lines
load("../county_tortoise_plotting_info.Robj")  # @gbradburd: how was this produced?

# and, coverages
coverages <- scan("../tortGen/exampleOutput/alleleCounts500kLoci.colmeans.txt")
dim(coverages) <- c(2,nrow(torts))
torts$coverage <- coverages[1,]
torts$mean.major.coverage <- coverages[2,]
torts$mean.major <- coverages[2,]/coverages[1,]


####
# ok now plot stuff

pdf(file="PCs-and-position.pdf", width=10, height=10, pointsize=10)
pairs( subset(torts,UTM_Zone!="12N")[c("PC1","PC2","PC3","PC4","Northing","Easting","coverage","mean.major")] )
dev.off()

# plots of PCs against each other, colored by northing and easting
pdf(file="PC-maps.pdf", width=10, height=8, pointsize=10)
    cols <- diverge_hcl(64)
    pairs( torts[c("PC1","PC2","PC3","PC4")], 
        lower.panel=function(x,y,...){points(x,y,bg=cols[cut(torts$Northing,breaks=64)], pch=21, cex=2, col=grey(.20) )}, 
        upper.panel=function(x,y,...){points(x,y,bg=cols[cut(torts$Easting,breaks=64)], pch=21, cex=2, col=grey(.20) )},
        main="Northing (above), Easting (below)" )
dev.off()

# maps, with PCs on
pdf(file="maps-with-PCs.pdf", width=10, height=8, pointsize=10)
    plot(elev,main="Elevation with tortoise IDs")
    lines(county_lines)
    text(tortoise_locations,labels=gsub("etort.","",torts$EM_Tort_ID))
    ncols <- 16
    cols <- adjustcolor(diverge_hcl(ncols),.7)
    for (pcs in c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")) {
        pcvec <- torts[[pcs]]
        pcfac <- cut( pcvec, breaks=ncols )
        plot(elev,main=pcs)
        lines(county_lines)
        points(tortoise_locations,pch=21,cex=2,col=grey(.2), bg=cols[pcfac])
        legend("topleft", legend=formatC(tapply(pcvec,pcfac,mean,na.rm=TRUE),digits=2), fill=cols)
    }
dev.off()

# check stuff for technical artifacts
if (interactive()) {
    for (x in grep("PC",names(torts),value=TRUE,invert=TRUE)) {
        layout(matrix(1:4,nrow=2))
        for (pc in c("PC1","PC2","PC3","PC4")) {
            plot( torts[[x]][torts$UTM_Zone!="12N"], torts[[pc]][torts$UTM_Zone!="12N"], main=paste(x,pc) )
        }
        if (is.null(locator(1))) break
    }
}


# look for association between PC1 and environmental variables
raster.covs <- cor( rasters, torts$PC1, use="pairwise" )
raster.covs <- as.vector( raster.covs ); names(raster.covs) <- colnames(rasters)

if (interactive()) {
    for (k in 1:ncol(rasters)) {
        plot( rasters[,k], torts$PC1, main=colnames(rasters)[k] )
        if (is.null(locator(1))) break
    }
}

# Does PC1 correlate with coverage?
pdf(file="PC1-and-coverage.pdf",width=6,height=6,pointsize=10)
    plot( PC1 ~ coverage, data=torts, type='n' )
    with(torts, text( coverage, PC1, labels=gsub("etort.","",EM_Tort_ID) ) )

    plot( PC1 ~ mean.major, data=torts, type='n' )
    with(torts, text( mean.major, PC1, labels=gsub("etort.","",EM_Tort_ID) ) )

    plot( PC1 ~ I(mean.major/coverage), data=torts, type='n' )
    with(torts, text( (mean.major/coverage), PC1, labels=gsub("etort.","",EM_Tort_ID) ) )
dev.off()


# how about correlations in coverage?
full.coverage.covmat <- as.matrix( read.table("../tortGen/exampleOutput/alleleCounts500k-raw-covmat.txt",header=FALSE) )
coverage.covmat <- full.coverage.covmat[ 2*(1:nrow(torts))-1, 2*(1:nrow(torts))-1 ]
minor.covmat <- full.coverage.covmat[ 2*(1:nrow(torts)), 2*(1:nrow(torts)) ]
coverage.minor.cormat <- cov2cor(full.coverage.covmat)[ 2*(1:nrow(torts))-1, 2*(1:nrow(torts)) ]

pdf(file="coverage-and-covariance.pdf",width=8,height=8,pointsize=10)
plot( as.vector(cov2cor(coverage.covmat)), as.vector(cov2cor(covmat)), col=adjustcolor(ifelse(row(covmat)==col(covmat),"red","black"),.5), pch=20 ,cex=.5, xlab="coverage correlation", ylab="allele frequency correlation" )
plot( as.vector(cov2cor(minor.covmat)), as.vector(cov2cor(covmat)), col=adjustcolor(ifelse(row(covmat)==col(covmat),"red","black"),.5), pch=20 ,cex=.5, xlab="major allele coverage correlation", ylab="allele frequency correlation" )
plot( as.vector(coverage.minor.cormat), as.vector(cov2cor(covmat)), col=adjustcolor(ifelse(row(covmat)==col(covmat),"red","black"),.5), pch=20 ,cex=.5, xlab="coverage-major allele cross-correlation", ylab="allele frequency correlation" )
dev.off()

if (FALSE) {

    ###
    # does not change after removing cluster in Ivanpah

    ivanpah <- ( torts$Location_ID %in% c("Ivanpah","ISEGS","Silver State") )

    sub.eig.covmat <- eigen( covmat[!ivanpah,!ivanpah] )
    rownames( sub.eig.covmat$vectors ) <- torts$EM_Tort_ID[ which(!ivanpah) ]

    for (k in 1:4) {
        torts[[paste("subPC",k,sep='')]] <- NA
        torts[[paste("subPC",k,sep='')]][which(!ivanpah)] <- sub.eig.covmat$vectors[,k]
    }

    pairs( subset(torts,UTM_Zone!="12N")[c("subPC1","subPC2","subPC3","subPC4","Northing","Easting")] )
}
