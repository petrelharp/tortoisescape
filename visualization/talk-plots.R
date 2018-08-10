#!/usr/bin/Rscript
require(methods)  # needed for get()
require(raster)
require(maps)
require(maptools)
require(rgdal)
require(rgeos)
require(TeachingDemos)

# elevation raster
dem <- raster("../visualization/dem_30.gri")

source("map-utils.R",chdir=TRUE)
shade <- get_shading(dem)
contours <- get_contours(dem)
states <- get_statelines(dem)
counties <- get_counties(dem)
ocean <- get_ocean(dem)

dist1.file <- "../tort_272_info/geog_distance.csv"
dist2.file <- "../tort_272_info/all_angsd_snps.pwp.csv"
indir <- "../tort_272_info"

dist1 <- read.csv(dist1.file,header=TRUE,stringsAsFactors=FALSE)
dist2 <- read.csv(dist2.file,header=TRUE,stringsAsFactors=FALSE)
# from Evan 1/31/15
dist2$pi <- dist2$raw_pi * ( 56575857 / 1898838430 )
# from Evan 2/10/15
dist2$years <- dist2$pi  / 2.064406e-8 / 2
# referenced in www.sciencedirect.com/science/article/pii/S0006320713003443
dist2$generations <- dist2$years / 25


# locations
tort.coord.obj <- load(file.path(indir,"geog_coords.RData"))
tort.coords <- get(tort.coord.obj)
tort.ids <- row.names(tort.coords)

# read in other info
pcs <- read.csv(file.path(indir,"pcs.csv"),header=TRUE,stringsAsFactors=FALSE)
stopifnot( all( tort.ids %in% pcs$etort ) )
base.cols <- c("blue","purple","red")  # north, south, between
pc.cols <- adjustcolor(base.cols,0.75)[ ifelse( pcs$PC1[match(tort.ids,pcs$etort)] > 0, 1, 2 ) ]

relatives <- ( (dist2$etort1==dist2$etort2) | ( dist2[,3] < quantile(subset(dist2,etort1==etort2)[,3],0.75) ) )

player <- function (main='') { 
    plot(dem,legend=FALSE,xlab="",ylab="",xaxt="n",yaxt="n",legend.mar=0,box=FALSE,main=main) 
    lines(counties,lwd=0.5)
}
pshade <- function (main='',samples=FALSE) { 
	plot(tort.coords,pch=NA,cex=0.7,ylab="",xlab="", main=main)
    plot( shade, col=adjustcolor(grey(seq(0,1,length.out=101)),0.25), legend=FALSE, add=TRUE )
    lines(counties,lwd=0.5)
    lines(states,lwd=2)
	scalebar(d=100000,xy=c(-1700000,y=-460000),below="meters",lwd=1.5,label="100km")
    if (samples) points(tort.coords,pch=20,cex=0.7)
}

png(file="everyone-pwp.png",width=5.5*144,height=5.5/1.7*144,pointsize=10,res=144)
    usethese <- !relatives
    north1 <- ( pcs$PC1[match(dist2$etort1[usethese],pcs$etort)] > 0 )
    north2 <- ( pcs$PC1[match(dist2$etort2[usethese],pcs$etort)] > 0 )
    layout(t(1:2))
    par(mar=c(2.5,2.5,0.5,0.5))
    player()
    points(tort.coords,pch=20,col=pc.cols)
    plot( dist1[usethese,3]/1000, dist2[usethese,3], pch=20, cex=.25, 
       col=adjustcolor("black",0.25), xlab="geog dist (km)", ylab="divergence",
       mgp=c(1.6,0.75,0) )
dev.off()

png(file="everyone-pwp-shaded.png",width=5.5*144,height=5.5/1.7*144,pointsize=10,res=144)
    usethese <- !relatives
    north1 <- ( pcs$PC1[match(dist2$etort1[usethese],pcs$etort)] > 0 )
    north2 <- ( pcs$PC1[match(dist2$etort2[usethese],pcs$etort)] > 0 )
    layout(t(1:2))
    par(mar=c(2.5,2.5,0.5,0.5))
    pshade()
    points(tort.coords,pch=20,col=pc.cols)
    plot( dist1[usethese,3]/1000, dist2[usethese,"years"], pch=20, cex=.25, 
       col=adjustcolor("black",0.25), xlab="geog dist (km)", ylab="divergence (years)",
       mgp=c(1.6,0.75,0) )
dev.off()

png(file="everyone-pwp-vertical.png",width=2.5*144,height=5*144,pointsize=10,res=144)
    usethese <- !relatives
    north1 <- ( pcs$PC1[match(dist2$etort1[usethese],pcs$etort)] > 0 )
    north2 <- ( pcs$PC1[match(dist2$etort2[usethese],pcs$etort)] > 0 )
    distcolors <- adjustcolor(base.cols,0.25)[ ifelse( north1&north2, 1, ifelse( (!north1)&(!north2), 2, 3 ) ) ]
    layout((1:2))
    par(mar=c(0,0,0.5,0))
    player()
    points(tort.coords,pch=20,col=pc.cols)
    par(mar=c(2.5,2.5,0.5,0.5))
    plot( dist1[usethese,3]/1000, dist2[usethese,3], pch=20, cex=.25, 
       col=distcolors, xlab="geog dist (km)", ylab="divergence",
       mgp=c(1.6,0.75,0) )
dev.off()

mindist <- min(dist2[!relatives,3])
sddist <- 3*sd(dist2[!relatives,3])
sfn <- function (x,max.cex=7) {
    max.cex/( 1 + exp( (x-mindist)/sddist ) )
}

for (tid in paste("etort-",c(285,240,35,273,57,229,191),sep='')) {
  png( file=paste("pwp_",tid,".png",sep=''), width=5.5*144, height=5.5/1.7*144, pointsize=10,res=144 )
    layout(t(1:2))
    par(mar=c(2.5,2.5,0.5,0.5))
    usethese <- ( dist2$etort1 != dist2$etort2 ) & ( ( dist2$etort1 == tid ) | ( dist2$etort2 == tid ) )
    otherone <- ifelse( dist2$etort1[usethese] == tid, dist2$etort2[usethese], dist2$etort1[usethese] )
    thiscolors <- pc.cols[ match(otherone,tort.ids) ]
    player()
    points(tort.coords[match(otherone,tort.ids)],pch=20,cex=sfn(dist2[,3][usethese]),col=thiscolors)
    points(tort.coords[match(tid,tort.ids)], pch="*", cex=4, col='black' )
    plot( dist1[!relatives,3]/1000, dist2[!relatives,3], pch=20, cex=.5, 
        col=adjustcolor("black",.25),
       xlab="geog dist (km)", ylab="divergence",
       mgp=c(1.6,0.75,0) )
    points( dist1[,3][usethese]/1000, dist2[,3][usethese], pch=20, col=thiscolors, cex=1.5 )
  dev.off()
}

# with shading
for (tid in paste("etort-",c(285,78,240,283,35,273,253,57,71,229,27,191),sep='')) {
  png( file=paste("pwp_",tid,"_shaded.png",sep=''), width=5.5*144, height=5.5/1.7*144, pointsize=10,res=144 )
    layout(t(1:2))
    par(mar=c(2.5,2.5,0.5,0.5)+.1)
    usethese <- ( dist2$etort1 != dist2$etort2 ) & ( ( dist2$etort1 == tid ) | ( dist2$etort2 == tid ) )
    otherone <- ifelse( dist2$etort1[usethese] == tid, dist2$etort2[usethese], dist2$etort1[usethese] )
    thiscolors <- pc.cols[ match(otherone,tort.ids) ]
    pshade()
    points(tort.coords[match(otherone,tort.ids)],pch=20,cex=sfn(dist2[,3][usethese]),col=thiscolors)
    points(tort.coords[match(tid,tort.ids)], pch="*", cex=4, col='black' )
    plot( dist1[!relatives,3]/1000, dist2[!relatives,"years"], pch=20, cex=.5, 
        col=adjustcolor("black",.25),
       xlab="geog dist (km)", ylab="divergence (years)",
       mgp=c(1.6,0.75,0) )
    points( dist1[,3][usethese]/1000, dist2[,"years"][usethese], pch=20, col=thiscolors, cex=1.5 )
  dev.off()
}

# for erik
tid <- "etort-285"
png(file="etort-285-simple.png", width=5.5*288, height=3*288, pointsize=10, res=288)
    layout(t(1:2))
    par(mar=c(2.5,2.5,1.0,0.5))
    usethese <- ( dist2$etort1 != dist2$etort2 ) & ( ( dist2$etort1 == tid ) | ( dist2$etort2 == tid ) ) & ( dist1$etort1 != "etort-1" ) & ( dist1$etort2 != "etort-1" )
    otherone <- ifelse( dist2$etort1[usethese] == tid, dist2$etort2[usethese], dist2$etort1[usethese] )
    badones <- (relatives | dist1$etort1 == "etort-1" | dist1$etort2 == "etort-1")
    plot(shade, col=adjustcolor(grey(seq(0,1,length.out=101)),0.5), legend=FALSE,
         xlim=c(-2000000,-1550000), ylim=c(-5e5, 0.5e5))
    plot(ocean, add = TRUE, col = "light blue")
    lines(contours,col=adjustcolor("black",0.2))
    lines(counties, lwd=0.5, col=adjustcolor("red",0.75))
    # player()
    points(tort.coords[match(otherone,tort.ids)],pch=20,cex=1, col=ifelse(otherone=="etort-1", NA, 'black'))
    points(tort.coords[match(tid,tort.ids)], pch=20, cex=2, col='blue' )
    plot( dist1[!badones,3]/1000, 1000 * dist2[!badones,3], pch=20, cex=.5, 
        col=adjustcolor("black",.25),
       xlab="geog dist (km)", ylab="divergence (per Kb)",
       mgp=c(1.6,0.75,0) )
    points( dist1[,3][usethese]/1000, 1000 * dist2[,3][usethese], pch=20, col="blue", cex=1.5 )
dev.off()

##########

source("../inference/resistance-fns.R")
require(raster)
require(rgdal)
show(load("../inference/habitat-only/nus_gt_three/setup.RData"))
show(load("../inference/habitat-only/nus_gt_three/inference-11254807_5.RData"))
load("../visualization/counties.Robj")  # provides counties



G@x <- update.G(trust.optim$argument[-1])
hts <- hitting.analytic(neighborhoods,G)

# see raster:::.rasterImagePlot
ph <- plot.ht.fn( layer.prefix="../geolayers/nussear/habitat-model/", nonmissing=nonmissing,sample.loc.file="../tort_272_info/geog_coords.RData")
clines <- spTransform( counties, CRSobj=CRS(proj4string(with(environment(ph),dem))))

png(file="example-hts.png",width=3*144,height=6*144,pointsize=10,res=144)
plot.inds <- c(16,98)
layout((1:2))
for (k in plot.inds) {
    ph( hts[,k], xaxt='n', yaxt='n', xlab='', ylab='', xlim=c(2e5,7.5e5), ylim=c(3.6e6,4.1e6), par.args=list(mar=c(0,0.1,0,0)+.5), zlim=c(0,20000), legend=FALSE, box=FALSE, legend.mar=0 )
    lines(clines)
}
dev.off()


# get tortoise coords and 
#	county lines in raster coordinate space
tort.coord.obj <- load("../tort_272_info/geog_coords.RData")
tort.coords <- get(tort.coord.obj)
raster_GCS_CRS_proj4 <- CRS(scan("../raster_GCS_CRS_proj4.txt",what="char",sep="\n"))

# load in elevation raster
dem <- raster("../visualization/dem_30.gri")

# define Ivanpah plotting region
x.min <- -1.825e6
x.max <- -1.672e6
y.min <- -3.0e05
y.max <- -1.5e05


# sample map on elevation
# png(file="sample_map_elev.png",res=200,width=10*200,height=5*200)
pdf(file="sample_map_elev.pdf",width=2,height=4,pointsize=10)
#quartz(width=10,height=5)
layout((1:2))
par(mar=c(0,0,1,0)+.5)
	plot(tort.coords,pch=20,cex=0.7,ylab="",xlab="")
    plot( shade, col=adjustcolor(grey(seq(0,1,length.out=101)),0.25), legend=FALSE, add=TRUE )
        points(tort.coords,pch=20,cex=0.7)
		lines(counties,lwd=0.5)
		lines(states,lwd=2)
		polygon(x=c(x.min,x.max,x.max,x.min,x.min),y=c(y.min,y.min,y.max,y.max,y.min),lty=2,lwd=0.5)
		box(lwd=3)
		mtext(side=3,font=2,text="Tortoise Sample Map",padj=-1)
	scalebar(d=100000,xy=c(-1650000,y=-460000),type="line",below="meters",lwd=1.5,label="100km")
	plot(tort.coords,pch=20,cex=0.7,ylab="",xlab="",xlim=c(x.min,x.max),ylim=c(y.min,y.max))
    plot( shade, col=adjustcolor(grey(seq(0,1,length.out=101)),0.5), legend=FALSE, add=TRUE )
        points(tort.coords,pch=20,cex=0.7)
		lines(counties,lwd=0.5)
		lines(states,lwd=2)
		box(lwd=3)
		mtext(side=3,font=2,text="Ivanpah Valley",padj=-1,line=0)
	scalebar(d=20000,xy=c(x.min+2000,y.max-15000),type="line",below="meters",lwd=1.1,label="20km")
dev.off()

# make a map of tortoise samples colored
#	by continuous position on PC1, 
#	on elevation map
tort.plotting.pch <- rep(20,nrow(tort.coords@coords))
tort.plotting.pch[
        ( pcs$PC2[match(tort.ids,pcs$etort)] > 0 )
        & ( pcs$PC1[match(tort.ids,pcs$etort)] < 0 )
    ] <- 17

#png(file="continuous_pc_colors_sample_map.png",res=288,width=2*288,height=3*288,pointsize=10)
pdf(file="continuous_pc_colors_sample_map.pdf",width=3.5,height=3,pointsize=10)
par(mar=c(0,0,2,0)+1)
pshade(main="PCs")
		points(tort.coords,pch=tort.plotting.pch,cex=0.7,col=pc.cols)
        subplot.coords <- cbind( x=grconvertX(.08+c(0,0.3),from='nfc',to='user'),
                                 y=grconvertY(.08+c(0,0.3),from='nfc',to='user') )
		rect(subplot.coords[1],
		     subplot.coords[3],
		     subplot.coords[2],
		     subplot.coords[4],col="white")
		subplot(fun = {		par(mgp=c(0.25,0,0));
                            plot(pcs$PC1,pcs$PC2,
								col=pc.cols,pch=tort.plotting.pch,cex=0.5,xaxt='n',yaxt='n',
                                xlab="PC1 (12%)",
                                ylab="PC2 (2%)",
                                ylim=c(-0.16,0.16)) ;
								abline(v=0,col="gray",lty=2) ;
								lines(x=c(-5,0),y=c(0,0),col="gray",lty=2)
						},
                        x=subplot.coords[,"x"], y=subplot.coords[,"y"] )
		# mtext(side=1,text="PC1 (6.9%)",cex=0.75,adj=0.215,padj=-3)
		# mtext(side=1,text="PC2 (1.4%)",cex=0.75,adj=-0.55,padj=-19.5,las=2)
		# mtext("meters",side=4,las=2,padj=-8.5,adj=-0.25)
	scalebar(d=100000,xy=c(-1570000,y=-560000),type="line",below="meters",lwd=1.5,label="100km")
dev.off()


