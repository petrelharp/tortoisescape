################################################################
################################################################
#	make figures for tortoise report
################################################################
################################################################

################################
#	pca maps
################################
require(raster)
require(maps)
require(maptools)
require(rgdal)
require(TeachingDemos)

if(file.exists("~/desktop/dropbox/tortoisescape")){
	setwd("~/desktop/dropbox/tortoisescape")
}
# get tortoise coords and 
#	county lines in raster coordinate space
load("tort_272_info/geog_coords.RData")
load("county_lines.Robj")
raster_GCS_CRS_proj4 <- CRS(scan("raster_GCS_CRS_proj4.txt",what="char",sep="\n"))

# and get state lines in bold
state.lines <- map(database="state",regions=c("California","Arizona","Nevada"),plot=FALSE)
state.lines.spobj <- map2SpatialLines(state.lines,proj4string=CRS("+proj=longlat"))
state.lines.spobj <- spTransform(state.lines.spobj,raster_GCS_CRS_proj4)

# load in elevation raster
if(file.exists("/Volumes/cooplab1/tortoises/geolayers/10x/crop_resampled_masked_aggregated_10x_dem_30.gri")){
	dem <- raster("/Volumes/cooplab1/tortoises/geolayers/10x/crop_resampled_masked_aggregated_10x_dem_30.gri")
}

# define Ivanpah plotting region
x.min <- -1.825e6
x.max <- -1.672e6
y.min <- -3.0e05
y.max <- -1.5e05

# make black & white sample map
png(file="tort_report_figures/sample.map.png",res=200,width=10*200,height=5*200)
#quartz(width=10,height=5)
par(mfrow=c(1,2))
	plot(tort.coords.rasterGCS,pch=20,cex=0.7,ylab="",xlab="")
		lines(county_lines,lwd=0.5)
		lines(state.lines.spobj,lwd=2)
		polygon(x=c(x.min,x.max,x.max,x.min,x.min),y=c(y.min,y.min,y.max,y.max,y.min),lty=2,lwd=0.5)
		box(lwd=3)
		mtext(side=3,font=2,text="Tortoise Sample Map",padj=-1)
	scalebar(d=100000,xy=c(-1590000,y=-460000),type="line",below="meters",lwd=1.5,label="100km")
	plot(tort.coords.rasterGCS,pch=20,cex=0.7,ylab="",xlab="",xlim=c(x.min,x.max),ylim=c(y.min,y.max))
		lines(county_lines,lwd=0.5)
		lines(state.lines.spobj,lwd=2)
		box(lwd=3)
		mtext(side=3,font=2,text="Ivanpah Sample Map",padj=-1)
	scalebar(d=20000,xy=c(x.min+2000,y.max-15000),type="line",below="meters",lwd=1.1,label="20km")
dev.off()

load("tort_272_info/geog_coords.RData")

png(file="tort_report_figures/black&white_272_sample.map.png",res=200,width=10*200,height=5*200)
#quartz(width=10,height=5)
par(mfrow=c(1,2))
	plot(tort.coords.272.rasterGCS,pch=20,cex=0.7,ylab="",xlab="",col="orange")
		points(tort.coords.rasterGCS,pch=20,cex=0.7)
		lines(county_lines,lwd=0.5)
		lines(state.lines.spobj,lwd=2)
		polygon(x=c(x.min,x.max,x.max,x.min,x.min),y=c(y.min,y.min,y.max,y.max,y.min),lty=2,lwd=0.5)
		box(lwd=3)
		mtext(side=3,font=2,text="Tortoise Sample Map",padj=-1)
	scalebar(d=100000,xy=c(-1590000,y=-460000),type="line",below="meters",lwd=1.5,label="100km")
	plot(tort.coords.272.rasterGCS,pch=20,cex=0.7,ylab="",xlab="",xlim=c(x.min,x.max),ylim=c(y.min,y.max),col="orange")
		points(tort.coords.rasterGCS,pch=20,cex=0.7)
		lines(county_lines,lwd=0.5)
		lines(state.lines.spobj,lwd=2)
		box(lwd=3)
		mtext(side=3,font=2,text="Ivanpah Sample Map",padj=-1)
	scalebar(d=20000,xy=c(x.min+2000,y.max-15000),type="line",below="meters",lwd=1.1,label="20km")
dev.off()


# do PCA on tortoise genetic covariance
#	and use it to make plotting color & symbol options
nind <- nrow(tort.coords.rasterGCS@coords)
covmat <- as.matrix(read.table("covmat/alleleCounts500kLoci-covmat.txt"))
pmat <- diag(nind) - 1/nind
eig.covmat <- eigen(pmat %*% covmat %*% pmat)
tort.plotting.colors.discrete <- rep("red",nrow(tort.coords.rasterGCS@coords))
	tort.plotting.colors.discrete[eig.covmat$vectors[,1] > 0] <- "purple"
tort.plotting.colors.continuous <- rev(rainbow(nind,start=4/6,end=6/6))[as.numeric(cut(eig.covmat$vectors[,1],nind))]
tort.plotting.pch <- rep(20,nrow(tort.coords.rasterGCS@coords))
	tort.plotting.pch[eig.covmat$vectors[,2] > 0 & eig.covmat$vectors[,1] < 0] <- 17
		xlab=,
		ylab=,


# make a map of tortoise samples colored
#	by N or S of Ivanpah on elevation raster map
png(file="tort_report_figures/discrete.pc.colors.sample.map.png",res=200,width=10*200,height=5*200)
#quartz(width=10,height=5)
par(mfrow=c(1,2),oma=c(1,1,1,1),mar=c(2,1,4,4))
#	plot(runif(10)) ; box(lwd=3)
plot(dem,ylab="",xlab="",main="Tortoise Sample Map",xaxt='n',yaxt='n')
		polygon(x=c(x.min,x.max,x.max,x.min,x.min),y=c(y.min,y.min,y.max,y.max,y.min),lty=2,lwd=0.5)
		#points(tort.coords.rasterGCS[match(c("etort-1","etort-16","etort-70"),row.names(tort.coords.rasterGCS)),],pch=1,col="green",cex=2)
			lines(county_lines,lwd=0.5)
			lines(state.lines.spobj,lwd=2)
		points(tort.coords.rasterGCS,pch=tort.plotting.pch,cex=0.7,col=tort.plotting.colors.discrete)
			box(lwd=3)
		rect(-2155000,-5.75e+05,-1872000,-4.2e+05,col="white")
		subplot(fun = {		plot(eig.covmat$vectors[,1],eig.covmat$vectors[,2],
								col=tort.plotting.colors.discrete,pch=tort.plotting.pch,cex=0.5,xaxt='n',yaxt='n',xlab="",ylab="") ; 
								abline(v=0,col="gray",lty=2) ;
								lines(x=c(-5,0),y=c(0,0),col="gray",lty=2)
						},
					x=c(-2155000,-1872000),y=c(-5.75e+05,-4.2e+05))
		mtext(side=1,text=paste("PC1 (",
					100*round(eig.covmat$values[1]/sum(eig.covmat$values),3)
					,"%)",sep=""),cex=0.75,adj=0.215,padj=-3)
		mtext(side=1,text=paste("PC2 (",
					100*round(eig.covmat$values[2]/sum(eig.covmat$values),3)
					,"%)",sep=""),cex=0.75,adj=-0.55,padj=-19.5,las=2)
		mtext("meters",side=4,las=2,padj=-8.5,adj=-0.25)
	scalebar(d=100000,xy=c(-1570000,y=-560000),type="line",below="meters",lwd=1.5,label="100km")
	zoom.ext <- extent(dem)
		zoom.ext@xmin <- x.min - 1e6
		zoom.ext@xmax <- x.max + 1e6
		zoom.ext@ymin <- y.min - 1.5e5
		zoom.ext@ymax <- y.max + 1.5e5
	plot(dem,ylab="",xlab="",main="Ivanpah Sample Map",ext=zoom.ext,xlim=c(x.min,x.max),ylim=c(y.min,y.max),asp=0.905,yaxt="n",xaxt="n")
		points(tort.coords.rasterGCS,pch=tort.plotting.pch,cex=0.7,col=tort.plotting.colors.discrete)
		#points(tort.coords.rasterGCS[match(c("etort-1","etort-16","etort-70"),row.names(tort.coords.rasterGCS)),],pch=1,col="green",cex=2)
		lines(county_lines,lwd=0.5)
		lines(state.lines.spobj,lwd=1)
		box(lwd=3)
	mtext("meters",side=4,las=2,padj=-8.5,adj=-0.25)
	scalebar(d=20000,xy=c(x.min+10000,y.max-15000),type="line",below="meters",lwd=1.1,label="20km")
dev.off()

# make a map of tortoise samples colored
#	by continuous position on PC1, 
#	on elevation raster map
png(file="tort_report_figures/continuous.pc.colors.sample.map.png",res=200,width=10*200,height=5*200)
#quartz(width=10,height=5)
par(mfrow=c(1,2),oma=c(1,1,1,1),mar=c(2,1,4,4))
#	plot(runif(10)) ; box(lwd=3)
plot(dem,ylab="",xlab="",main="Tortoise Sample Map",xaxt='n',yaxt='n')
		polygon(x=c(x.min,x.max,x.max,x.min,x.min),y=c(y.min,y.min,y.max,y.max,y.min),lty=2,lwd=0.5)
		#points(tort.coords.rasterGCS[match(c("etort-1","etort-16","etort-70"),row.names(tort.coords.rasterGCS)),],pch=1,col="green",cex=2)
			lines(county_lines,lwd=0.5)
			lines(state.lines.spobj,lwd=2)
		points(tort.coords.rasterGCS,pch=tort.plotting.pch,cex=0.7,col=tort.plotting.colors.continuous)
			box(lwd=3)
		rect(-2155000,-5.75e+05,-1872000,-4.2e+05,col="white")
		subplot(fun = {		plot(eig.covmat$vectors[,1],eig.covmat$vectors[,2],
								col=tort.plotting.colors.continuous,pch=tort.plotting.pch,cex=0.5,xaxt='n',yaxt='n',xlab="",ylab="",ylim=c(-0.16,0.16)) ;
								abline(v=0,col="gray",lty=2) ;
								lines(x=c(-5,0),y=c(0,0),col="gray",lty=2)
						},
					x=c(-2155000,-1872000),y=c(-5.75e+05,-4.2e+05))
		mtext(side=1,text=paste("PC1 (",
					100*round(eig.covmat$values[1]/sum(eig.covmat$values),3)
					,"%)",sep=""),cex=0.75,adj=0.215,padj=-3)
		mtext(side=1,text=paste("PC2 (",
					100*round(eig.covmat$values[2]/sum(eig.covmat$values),3)
					,"%)",sep=""),cex=0.75,adj=-0.55,padj=-19.5,las=2)
		mtext("meters",side=4,las=2,padj=-8.5,adj=-0.25)
	scalebar(d=100000,xy=c(-1570000,y=-560000),type="line",below="meters",lwd=1.5,label="100km")
	zoom.ext <- extent(dem)
		zoom.ext@xmin <- x.min - 1e6
		zoom.ext@xmax <- x.max + 1e6
		zoom.ext@ymin <- y.min - 1.5e5
		zoom.ext@ymax <- y.max + 1.5e5
	plot(dem,ylab="",xlab="",main="Ivanpah Sample Map",ext=zoom.ext,xlim=c(x.min,x.max),ylim=c(y.min,y.max),asp=0.905,yaxt="n",xaxt="n")
		points(tort.coords.rasterGCS,pch=tort.plotting.pch,cex=0.7,col=tort.plotting.colors.continuous)
		#points(tort.coords.rasterGCS[match(c("etort-1","etort-16","etort-70"),row.names(tort.coords.rasterGCS)),],pch=1,col="green",cex=2)
		lines(county_lines,lwd=0.5)
		lines(state.lines.spobj,lwd=1)
		box(lwd=3)
	mtext("meters",side=4,las=2,padj=-8.5,adj=-0.25)
	scalebar(d=20000,xy=c(x.min+10000,y.max-15000),type="line",below="meters",lwd=1.1,label="20km")
dev.off()

png(file="tort_report_figures/PC1_map.png",res=200,width=5*200,height=5*200)
	plot(eig.covmat$vectors[,1],
		eig.covmat$vectors[,2],
		col=tort.plotting.colors.continuous,
		pch=tort.plotting.pch,
		cex=1.5,main="Tortoise Principal Components Analysis",
		xlab=paste("PC1 (",
					100*round(eig.covmat$values[1]/sum(eig.covmat$values),3)
					,"%)",sep=""),
		ylab=paste("PC2 (",
					100*round(eig.covmat$values[2]/sum(eig.covmat$values),3)
					,"%)",sep=""),
		cex.axis=0.8)
dev.off()

png(file="tort_report_figures/PC1_map_partitioned.png",res=200,width=5*200,height=5*200)
	#quartz(width=5,height=5)
	plot(eig.covmat$vectors[,1],
		eig.covmat$vectors[,2],
		col=tort.plotting.colors.continuous,
		pch=tort.plotting.pch,
		cex=1.5,main="Tortoise Principal Components Analysis",
		xlab=paste("PC1 (",
					100*round(eig.covmat$values[1]/sum(eig.covmat$values),3)
					,"%)",sep=""),
		ylab=paste("PC2 (",
					100*round(eig.covmat$values[2]/sum(eig.covmat$values),3)
					,"%)",sep=""),
		cex.axis=0.8,
		ylim=c(-0.16,0.16))
		abline(v=0,col="gray",lty=2)
		lines(x=c(-5,0),y=c(0,0),col="gray",lty=2)
dev.off()


if(FALSE){
quartz(width=10,height=5)
layout(matrix(c(1,2,1,1,3,3,3,3), 2, 4))
#	plot(runif(10)) ; box(lwd=3)
plot(dem,ylab="",xlab="",main="Map of Tortoise Samples",xaxt='n',yaxt='n')
		points(tort.coords.rasterGCS,pch=20,cex=0.7,col=tort.plotting.colors)
			lines(county_lines,lwd=0.5)
			lines(state.lines.spobj,lwd=2)
			box(lwd=3)
	plot(eig.covmat$vectors[,1],eig.covmat$vectors[,2],
			col=tort.plotting.colors,pch=20,cex=0.5,xaxt='n',yaxt='n')
		box(lwd=1)
	zoom.ext <- extent(dem)
		zoom.ext@xmin <- x.min - 1e6
		zoom.ext@xmax <- x.max + 1e6
		zoom.ext@ymin <- y.min - 1e5
		zoom.ext@ymax <- y.max + 1e5
	plot(dem,ylab="",xlab="",main="Map of Tortoise Samples",ext=zoom.ext,xlim=c(x.min,x.max),ylim=c(y.min,y.max),asp=0.89,yaxt="",xaxt="")
		points(tort.coords.rasterGCS,pch=20,cex=0.7,col=tort.plotting.colors)
		lines(county_lines,lwd=0.5)
		lines(state.lines.spobj,lwd=1)
		box(lwd=3)

quartz(width=10,height=5)
par(mfrow=c(1,2))
#	plot(runif(10)) ; box(lwd=3)
plot(dem,ylab="",xlab="",main="Map of Tortoise Samples",xaxt='n',yaxt='n')
		points(tort.coords.rasterGCS,pch=20,cex=0.7,col=tort.plotting.colors)
			lines(county_lines,lwd=0.5)
			lines(state.lines.spobj,lwd=2)
			box(lwd=3)
		polygon(x=c(),y=c(),
		par(fig=c(0.06,0.25,0.18,0.45),new=TRUE)
			plot.window(bg="white")
			plot(eig.covmat$vectors[,1],eig.covmat$vectors[,2],
					col=tort.plotting.colors,pch=20,cex=0.5,xaxt='n',yaxt='n')
			box(lwd=1)
	zoom.ext <- extent(dem)
		zoom.ext@xmin <- x.min - 1e6
		zoom.ext@xmax <- x.max + 1e6
		zoom.ext@ymin <- y.min - 1e5
		zoom.ext@ymax <- y.max + 1e5
	plot(dem,ylab="",xlab="",main="Map of Tortoise Samples",ext=zoom.ext,xlim=c(x.min,x.max),ylim=c(y.min,y.max),asp=0.89,yaxt="",xaxt="")
		points(tort.coords.rasterGCS,pch=20,cex=0.7,col=tort.plotting.colors)
		lines(county_lines,lwd=0.5)
		lines(state.lines.spobj,lwd=1)
		box(lwd=3)
}
