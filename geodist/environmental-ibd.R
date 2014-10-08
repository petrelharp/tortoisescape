#!/usr/bin/Rscript


if (!interactive()) {
    subdir <- commandArgs(TRUE)[1]
    layer.file <- commandArgs(TRUE)[2]
} else {
    subdir <- "10x"
    layer.file <- "../inference/twentyfour-raster-list"
}
layer.names <- scan(layer.file,what="char")

require(raster)
require(colorspace)
require(MASS)

# geographic distance
torts <- read.csv("../1st_180_torts.csv",header=TRUE,stringsAsFactors=FALSE)
nind <- nrow(torts)
tort.dist.table <- read.table("../1st180_pairwise_distances_sorted_redundancy_removed.txt",header=TRUE)
tort.dists <- numeric(nind^2); dim(tort.dists) <- c(nind,nind)
tort.dists[ cbind( match(tort.dist.table$etort1,torts$EM_Tort_ID), match(tort.dist.table$etort2,torts$EM_Tort_ID) ) ] <- tort.dist.table$DISTANCE
tort.dists <- tort.dists + t(tort.dists)

# environmental distance
edist.list <- lapply( layer.names, function (layer.name) {
        x <- read.table(paste(subdir,"/",layer.name,"-mean-edist.tsv",sep=''),sep='\t',header=TRUE)
        colnames(x)[3] <- layer.name
        x
    } )
edists <- edist.list[[1]]
for (k in setdiff( seq_along(edist.list), 1 ) ) {
    edists <- merge( edists, edist.list[[k]], by=c("tort1","tort2") )
}

dists <- read.csv("../pairwise-normalized-pi.csv") # has DISTANCE, pi, and npi
dists <- merge( dists, edists, by.x=c("etort1","etort2"), by.y=c("tort1","tort2") )
dists$etort1 <- factor( dists$etort1 , levels=torts$EM_Tort_ID )
dists$etort2 <- factor( dists$etort2 , levels=torts$EM_Tort_ID )

##
# look at nearby distances
dist.cutoff <- 100
dists <- subset(dists, DISTANCE<dist.cutoff)

edist.lms <- lapply( layer.names, function (layer.name) {
                z <- dists[[ layer.name ]]
                lm( dists$pi ~ z )
    } )
names(edist.lms) <- layer.names

edist.coefs <- t( rbind( sapply( edist.lms, coef ), 
        r.squared=sapply( lapply( edist.lms, summary ), "[[", "r.squared" ),
        p.value=sapply( lapply( lapply( edist.lms, anova ), "[[", "Pr(>F)" ), "[", 1 ) 
    ) )
edist.coefs <- edist.coefs[ order(edist.coefs[,3]), ]

write.table( edist.coefs, file="envdist-coefficients.tsv", sep='\t', quote=FALSE )

pdf(file="envdist-correlations.pdf", width=12, height=8, pointsize=10 )
layout( matrix(1:24,nrow=4) )
par(mar=c(0,0,2,0)+.1)
for (k in seq_along(layer.names)) {
    z <- dists[[ layer.names[k] ]]
    plot( z, dists$pi, main=layer.names[k], xlab='', ylab='', xaxt='n', yaxt='n', pch=20, cex=.5, 
        col=adjustcolor("black",.5), xlim=quantile(z,c(0,.98),na.rm=TRUE) )
    abline(coef(edist.lms[[k]]),col='red',lwd=2)
}
dev.off()

if (FALSE) {

pcs <- read.csv("../covmat/tort-PCs.csv",stringsAsFactors=FALSE)
names(pcs)[1] <- "etort"
pcs$etort <- factor( pcs$etort, levels=torts$EM_Tort_ID )
stopifnot( all(pcs$etort==torts$EM_Tort_ID) )

ncols <- 16
cols <- adjustcolor( rev(sequential_hcl(ncols,h=40)), .5 )
efac <- cut(dists$mean_slope_30,ncols)

with(pcs, plot( PC1, PC2 ) )
for (lev in 1:nlevels(efac)) {
    with( subset( dists, DISTANCE < 40 & efac==levels(efac)[lev]), segments( 
            x0=pcs$PC1[etort1], y0=pcs$PC2[etort1], 
            x1=pcs$PC1[etort2], y1=pcs$PC2[etort2],
            col=cols[lev]
            ) )
}
legend("bottomright", fill=cols[2*(1:(ncols/2))], legend=formatC(tapply(dists$mean_slope_30,cut(dists$mean_slope_30,ncols),mean),digits=2)[2*(1:(ncols/2))])


with( subset(dists,pcs$PC1[match(etort1,pcs$etort)]*pcs$PC1[match(etort2,pcs$etort)]>=0),
    hist( mean_slope_30, col=adjustcolor("red",.5) , xlim=range(dists$mean_slope_30),freq=FALSE )
    )
with( subset(dists,pcs$PC1[match(etort1,pcs$etort)]*pcs$PC1[match(etort2,pcs$etort)]<0),
    hist( mean_slope_30, col=adjustcolor("blue",.5), add=TRUE,freq=FALSE )
    )


layer <- raster("../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_dem_30")
load("../tort.coords.rasterGCS.Robj")
coord.names <- rownames(tort.coords.rasterGCS@coords)

plot(layer,main="adjustment factor")
text(tort.coords.rasterGCS, gsub("etort.","",coord.names),cex=.5)

}
