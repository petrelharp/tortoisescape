#!/usr/bin/Rscript
setwd("..")  # move to base of tortoisescape

load("tort_180_info/tort.coords.rasterGCS.Robj")

# converting these pi values to divergence:
# Evan: "So the total size of the first hundred scaffolds is 285,214,272bp. In those 285,214,272bp we found a total of 5,319,387 SNPs with that angsd run."
pi.factor <- 5319387/285214272

require(raster)
layer <- raster("geolayers/TIFF/10x/crop_resampled_masked_aggregated_10x_dem_30.gri")

torts <- read.csv("1st_180_torts.csv",header=TRUE,stringsAsFactors=FALSE)
torts$EM_Tort_ID <- factor( torts$EM_Tort_ID , levels=torts$EM_Tort_ID )
nind <- nrow(torts)

dists <- read.csv("pairwise-normalized-pi.csv",header=TRUE,stringsAsFactors=FALSE) # had DISTANCE, pi, and npi
dists$etort1 <- factor( dists$etort1 , levels=torts$EM_Tort_ID )
dists$etort2 <- factor( dists$etort2 , levels=torts$EM_Tort_ID )
dists$pi <- dists$pi * pi.factor
dists$npi <- dists$npi * pi.factor

pcs <- read.csv("covmat/tort-PCs.csv",stringsAsFactors=FALSE)
pcs$X <- factor( pcs$X, levels=torts$EM_Tort_ID )
pc.cols <- adjustcolor( ifelse( pcs$PC1 > 0, "blue", "purple" ), .75 )
stopifnot( all(pcs$X==torts$EM_Tort_ID) )

# overall
png( file="pngs/everyone-ibd.png", width=12*144, height=4*144, pointsize=10, res=144 )
    thiscolors <- with(dists, ifelse( ( pcs$PC1[match(etort1,pcs$X)] > 0 ), 
                    ifelse ( ( pcs$PC1[match(etort2,pcs$X)] > 0 ), "green", "black" ),
                    ifelse ( ( pcs$PC1[match(etort2,pcs$X)] > 0 ), "black", "blue" ) ) )
    layout(t(1:3))
    plot(layer)
    points(tort.coords.rasterGCS,pch=20,cex=1,col=pc.cols)
    plot( dists$DISTANCE, dists$pi, pch=20, cex=.5, 
        col=adjustcolor("black",.25), xlab="geographic distance (km)", ylab="raw pairwise divergence" )
    points( dists$DISTANCE, dists$pi, pch=20, col=thiscolors, cex=1.5 )
    plot( dists$DISTANCE, dists$npi, pch=20, cex=.5, 
        col=adjustcolor("black",.25), xlab="geographic distance (km)", ylab="adjusted pairwise divergence" )
    points( dists$DISTANCE, dists$npi, pch=20, col=thiscolors, cex=1.5 )
dev.off()

# highlighted by individual
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


# pcs
pcs <- read.csv("covmat/tort-PCs.csv",stringsAsFactors=FALSE)
names(pcs)[1] <- "etort"
pcs$etort <- factor( pcs$etort, levels=torts$EM_Tort_ID )
stopifnot( all(pcs$etort==torts$EM_Tort_ID) )
torts <- cbind(torts,pcs[,-1])

# same-different groups
dists$pcgroup <- with(dists, ( (torts$PC1[match(etort1,torts$EM_Tort_ID)]*torts$PC1[match(etort2,torts$EM_Tort_ID)]) > 0 ) )
dists$PC1 <- with(dists, torts$PC1[match(etort1,torts$EM_Tort_ID)] - torts$PC1[match(etort2,torts$EM_Tort_ID)] )


##
# only within Ivanpah

torts$Location_ID <- factor(torts$Location_ID)
ivanpah_locs <- c("Silver State", "ISEGS", "Ivanpah")
ivanpah_torts <- torts$EM_Tort_ID[ torts$Location_ID %in% ivanpah_locs ]
ivanpah_coords <- tort.coords.rasterGCS[ torts$Location_ID %in% ivanpah_locs ]

ivanpah_loc_factor <- levels(torts$Location_ID)[ torts$Location_ID[ torts$Location_ID %in% ivanpah_locs ] ]
ivanpah_loc_factor <- factor( paste( ivanpah_loc_factor, ifelse( ivanpah_loc_factor=="ISEGS", 
                    ifelse( torts$Easting[torts$Location_ID %in% ivanpah_locs] < 645500, "_east", "_west" ),
                    "" ), sep='' ) )

ivanpah_dists <- subset( dists, ( etort1 %in% ivanpah_torts ) & ( etort2 %in% ivanpah_torts ) )
ivanpah_dists$loc1 <- ivanpah_loc_factor[ match(ivanpah_dists$etort1,ivanpah_torts) ]
ivanpah_dists$loc2 <- ivanpah_loc_factor[ match(ivanpah_dists$etort2,ivanpah_torts) ]
ivanpah_dists$col <- with(ivanpah_dists, ifelse( loc1==loc2, "green", "black") )

ivanpah_lm <- lm( pi ~ DISTANCE, data=ivanpah_dists )
ivanpah_within_lm <- lm( pi ~ DISTANCE, data=subset(ivanpah_dists,loc1==loc2) )

infl <- function (x,fac=.5) { as.vector(mean(x) + (1+fac)*(x-mean(x))) }

require(plotrix)

pdf(file="visualization/plots/ivanpah-ibd.pdf", width=10, height=6, pointsize=10)
layout(t(1:2))
plot(layer, xlim=infl(ivanpah_coords@bbox["coords.x1",],.8), ylim=infl(ivanpah_coords@bbox["coords.x2",]), main="Ivanpah valley", xaxt='n', yaxt='n' )
points( ivanpah_coords, pch=20 ) #, col=torts$Location_ID[torts$Location_ID %in% ivanpah_locs] )
with(ivanpah_dists, plot( DISTANCE, pi, xlab="geographic distance (km)", ylab="genetic distance (divergence)", pch=20, main=paste("IBD for", length(ivanpah_torts), "tortoises in Ivanpah") ) ) # , col=col ) )
abline( coef( ivanpah_lm ), col='red' )
# abline( coef( ivanpah_within_lm ), col='green' )
addtable2plot( "bottomright", table=round(summary(ivanpah_lm)$coefficients,digits=6), display.rownames=TRUE)
plot(layer, xlim=infl(ivanpah_coords@bbox["coords.x1",],.8), ylim=infl(ivanpah_coords@bbox["coords.x2",]), main="Ivanpah valley", xaxt='n', yaxt='n' )
points( ivanpah_coords, pch=20 , col=ivanpah_loc_factor )
with(ivanpah_dists, plot( DISTANCE, pi, xlab="geographic distance (km)", ylab="genetic distance (divergence)", pch=20, main=paste("IBD for", length(ivanpah_torts), "tortoises in Ivanpah"), col=col ) )
abline( coef( ivanpah_within_lm ), col='green' )
addtable2plot( "bottomright", table=round(summary(ivanpah_within_lm)$coefficients,digits=6), display.rownames=TRUE)
dev.off()


png(file="tort_report_figures/within-ivanpah-ibd.png",res=200,width=6*200,height=4*200)
layout(t(1:2))
    opar <- par(mar=c(5,4,3,1)+.1)
    with(ivanpah_dists, plot( DISTANCE, pi, 
			pch=20, cex=0.3, 
            xlab="geographic distance (km)", ylab="genetic distance (divergence)") )
    title( main="Isolation by Distance", line=1 )
    abline( coef( ivanpah_lm ), col='red' )
    par(mar=c(5,0,3,1)+.1)
    plot(layer, xlim=infl(ivanpah_coords@bbox["coords.x1",],1.2), ylim=infl(ivanpah_coords@bbox["coords.x2",]), xaxt='n', yaxt='n', 
        smallplot=c(0.68,0.71,0.2845,0.8155), 
        bigplot=c(0,0.64,0.255,0.845), 
        legend.lab="elevation (m)" )
    title( main="in the Ivanpah valley", line=1 )
    points( ivanpah_coords, pch=20 ) #, col=torts$Location_ID[torts$Location_ID %in% ivanpah_locs] )
    scalebar(1e4,label=c("","10km",""))
dev.off()
