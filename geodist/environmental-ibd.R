#!/usr/bin/Rscript


# if (!interactive()) {
#     subdir <- commandArgs(TRUE)[1]
#     layer.file <- commandArgs(TRUE)[2]
# } else {
    subdir <- "10x"
    layer.file <- "../inference/twentyfour-raster-list"
# }
layer.names <- scan(layer.file,what="char")

require(raster)
require(colorspace)
require(MASS)

# converting these pi values to divergence:
# Evan: "So the total size of the first hundred scaffolds is 285,214,272bp. In those 285,214,272bp we found a total of 5,319,387 SNPs with that angsd run."
pi.factor <- 5319387/285214272

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
dists$pi <- dists$pi * pi.factor
dists$npi <- dists$npi * pi.factor
dists <- merge( dists, edists, by.x=c("etort1","etort2"), by.y=c("tort1","tort2") )
dists$etort1 <- factor( dists$etort1 , levels=torts$EM_Tort_ID )
dists$etort2 <- factor( dists$etort2 , levels=torts$EM_Tort_ID )

# pcs
pcs <- read.csv("../covmat/tort-PCs.csv",stringsAsFactors=FALSE)
names(pcs)[1] <- "etort"
pcs$etort <- factor( pcs$etort, levels=torts$EM_Tort_ID )
stopifnot( all(pcs$etort==torts$EM_Tort_ID) )
torts <- cbind(torts,pcs[,-1])

# same-different groups
dists$pcgroup <- with(dists, ( (torts$PC1[match(etort1,torts$EM_Tort_ID)]*torts$PC1[match(etort2,torts$EM_Tort_ID)]) > 0 ) )
dists$PC1 <- with(dists, torts$PC1[match(etort1,torts$EM_Tort_ID)] - torts$PC1[match(etort2,torts$EM_Tort_ID)] )

# residual pi after accounting for geographical distance and PC1 distance
dists$resid.pi <- with(dists, resid( lm( pi ~ DISTANCE + PC1 ) ) )
##
# look at nearby distances
dist.cutoff <- 100
usethese <- with( dists, DISTANCE<dist.cutoff & pcgroup )

edist.lms <- lapply( layer.names, function (layer.name) {
                z <- subset(dists,usethese)[[ layer.name ]]
                lm( subset(dists,usethese)$pi ~ z )
    } )
names(edist.lms) <- layer.names

edist.coefs <- t( rbind( sapply( edist.lms, coef ), 
        r.squared=sapply( lapply( edist.lms, summary ), "[[", "r.squared" ),
        p.value=sapply( lapply( lapply( edist.lms, anova ), "[[", "Pr(>F)" ), "[", 1 ) 
    ) )
edist.coefs <- edist.coefs[ order(edist.coefs[,3]), ]

edist.resid.lms <- lapply( layer.names, function (layer.name) {
                z <- subset(dists,usethese)[[ layer.name ]]
                lm( subset(dists,usethese)$resid.pi ~ z )
    } )
names(edist.resid.lms) <- layer.names

edist.resid.coefs <- t( rbind( sapply( edist.resid.lms, coef ), 
        r.squared=sapply( lapply( edist.resid.lms, summary ), "[[", "r.squared" ),
        p.value=sapply( lapply( lapply( edist.resid.lms, anova ), "[[", "Pr(>F)" ), "[", 1 ) 
    ) )
layer.order <- rev(order(edist.resid.coefs[,3]))
edist.resid.coefs <- edist.resid.coefs[ layer.order, ]

options(scipen=5)
tmp.tex <- apply(edist.resid.coefs,2,format,digits=2)
write.table( tmp.text, file="envdist-coefficients.tsv", sep='\t', quote=FALSE )

# pdf(file="envdist-correlations.pdf", width=12, height=8, pointsize=10 )
png(file="envdist-correlations.png", width=12*288, height=8*288, res=288, pointsize=10 )
layout( matrix(1:24,nrow=4,byrow=TRUE) )
par(mar=c(0,0,2,0)+.1)
# cols <- adjustcolor("black",0.5)
cols <- adjustcolor( ifelse( abs(subset(dists,usethese)$PC1)<0.1, "black", "red" ), 0.5 )
# for (k in seq_along(layer.names)) {
#     z <- subset(dists,usethese)[[ layer.names[k] ]]
#     plot( z, subset(dists,usethese)$pi, main=layer.names[k], xlab='', ylab='', xaxt='n', yaxt='n', pch=20, cex=.5, 
#         col=cols, xlim=quantile(z,c(0,.98),na.rm=TRUE) )
#     abline(coef(edist.lms[[k]]),col='red',lwd=2)
# }
for (k in seq_along(layer.names)[layer.order]) {
    z <- subset(dists,usethese)[[ layer.names[k] ]]
    plot( z, subset(dists,usethese)$resid.pi, main=layer.names[k], xlab='', ylab='', xaxt='n', yaxt='n', pch=20, cex=.5, 
        col=cols, xlim=quantile(z,c(0,.98),na.rm=TRUE) )
    abline(coef(edist.resid.lms[[k]]),col='red',lwd=2)
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
