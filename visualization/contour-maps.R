#!/usr/bin/Rscript
library(methods)  # required for get( ) below... ???

usage <- "Makes many images of contour maps for a single data frame of pairwise comparisons: one for each thing being compared.
Usage:
    Rscript distance-maps.R (genetic distance file) (sample info directory) (output directory)
"

argvec <- if (interactive()) { scan(what="char") } else { commandArgs(TRUE) }

if (length(argvec)<3) { stop(usage) }

dist.file <- argvec[1]
indir <- argvec[2]
outdir <- argvec[3]

# dist.file <- "all_angsd_snps.pwp.csv"
dists <- read.csv(dist.file,header=TRUE,stringsAsFactors=FALSE)
dir.create(outdir,showWarnings=FALSE)

# locations
coord.obj <- load(file.path(indir,"geog_coords.RData"))
coords <- get(coord.obj)
tort.ids <- row.names(coords)
# read in other info
pcs <- read.csv(file.path(indir,"pcs.csv"),header=TRUE,stringsAsFactors=FALSE)
stopifnot( all( tort.ids %in% pcs$etort ) )
pc.cols <- adjustcolor( ifelse( pcs$PC1[match(tort.ids,pcs$etort)] > 0, "blue", "purple" ), .75 )
require(raster)
layer <- raster("../visualization/dem_30")
player <- function (main='',...) { plot(layer,legend=FALSE,xlab="",ylab="",xaxt="n",yaxt="n",legend.mar=0,box=FALSE,main=main,...) }

require(TSP)
xy <- coordinates(coords)
etsp <- ETSP( xy, labels=rownames(xy) )
tour <- solve_TSP( etsp, method="linkern" )
tour.labels <- t(outer(letters,letters,paste,sep=''))[seq_len(length(tour))]

library(fields)

xy <- coords@coords
dists$x1 <- xy[match(dists$etort1,rownames(xy)),1]
dists$y1 <- xy[match(dists$etort1,rownames(xy)),2]
dists$x2 <- xy[match(dists$etort2,rownames(xy)),1]
dists$y2 <- xy[match(dists$etort2,rownames(xy)),2]

fatten <- function (x,fac,mu=median(x,na.rm=TRUE)) { mu + (x-mu)*(1+fac) }
other <- function (etort) {
    with( subset(dists,(etort1==etort | etort2==etort) & (etort1!=etort2)),
            data.frame( 
                       etort=ifelse(etort1==etort,etort2,etort1),
                       pi=pi,
                       x=ifelse(etort1==etort,x2,x1),
                       y=ifelse(etort1==etort,y2,y1)
         ) )
}


xvals <- seq(min(xy[,1]),max(xy[,1]),length.out=100)
yvals <- seq(min(xy[,2]),max(xy[,2]),length.out=80)
xygrid <- expand.grid( x=xvals, y=yvals )
zlevels <- with(subset(dists,etort1!=etort2),seq(0.95*min(pi),1.05*max(pi),length.out=60))

lopred <- function (etort,...) {
    # loess
    xyl <- loess( pi ~ x*y,data=other(etort), span=0.1 )
    predict( xyl, newdata=xygrid )
}

krpred <- function (etort,theta=3e4,...) {
    # Kriging with covariance width theta
    xyz <- other(etort)
    xyk <- fields::Krig(x=xyz[,c("x","y")],Y=xyz$pi,theta=theta,...)
    matrix( predict( xyk, xygrid ), nrow=length(xvals), ncol=length(yvals) )
}

pi.contour <- function (etort,predfn=lopred,...) {
    pred <- predfn(etort,...)
    player(main=etort,xlim=range(xygrid$x),ylim=range(xygrid$y))
    points(coords,pch=20,cex=0.5)
    contour( xvals, yvals, z=pred, add=TRUE, levels=zlevels, nlevels=40, col=adjustcolor("black",0.5) )
    points( coords[etort], cex=3, col=adjustcolor('red',0.75), pch=20 )
}

# pi.contour("etort-42", predfn=krpred )

for (k in seq_along(tort.ids)) {
    tid <- tort.ids[tour[k]]
    cat(tid,"\n")
  png( file=file.path(outdir,paste(tour.labels[k], "_",gsub("[^0-9a-z-]","_",tid),".png",sep='')), width=4*144, height=4*144, pointsize=10, res=144 )
    pi.contour( tid, predfn=krpred )
  dev.off()
}
