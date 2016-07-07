
library(rgdal)
library(raster)
library(sp)
library(rgeos)
library(spdep)
library(colorspace)
library(maps)
source("../../visualization/map-utils.R",chdir=TRUE)
source("watersheds-fns.R")

# stuff for plotting
tort.coords.obj <- load("../../tort_272_info/geog_coords.RData")
tort.coords <- get(tort.coords.obj)


dem <- get_elev()
conts <- get_contours(dem)
shade <- get_shading(dem)

wbd.list <- lapply( 2*(1:5), function (wbd_num) {
            gClip( 
                  spTransform( 
                               readOGR("WBD/WBD_National_931.gdb",paste0("WBDHU",wbd_num)), 
                           CRS(proj4string(dem)) ), 
                     dem )
    } )


hierarchies <- c( list( rep(1,length(wbd.list[[1]])) ),
                  lapply( seq_along(wbd.list)[-1], function (k) {
                      hierarchy( wbd.list[[k]], wbd.list[[k-1]] ) } ) 
                 )


cols <- color_hierarchy( hierarchies )

plot_wbd <- function (k, labels=seq_along(wbd.list[[k]]), do.torts=TRUE,...) {
    plot( wbd.list[[k]], col=cols[[k]], ...,
             xlim=.expand(extent(tort.coords)[1:2],.1), 
             ylim=.expand(extent(tort.coords)[3:4],.1) )
    plot(shade, col=adjustcolor(grey(seq(0,1,length.out=101)),0.5), legend=FALSE, add=TRUE )
    lines(conts,col=adjustcolor("black",0.25))
    if (do.torts) points(tort.coords, pch=20, cex=1)
    xy <- gCentroid( wbd.list[[k]], byid=TRUE )
    text(xy,labels)
}

# check for contiguity
layout(matrix(1:8,nrow=2))
for (k in 2:5) {
    plot_wbd(k-1, col=cols[[k-1]], do.torts=FALSE)
    plot_wbd(k, labels=hierarchies[[k]], col=cols[[k]], do.torts=FALSE)
}


layout(1)

k <- 4
    plot_wbd(k, col=cols[[k]], labels=NULL, 
             xlim=.expand(extent(tort.coords)[1:2],.1), 
             ylim=.expand(extent(tort.coords)[3:4],.1) )
    plot(shade, col=adjustcolor(grey(seq(0,1,length.out=101)),0.5), legend=FALSE, add=TRUE )
    lines(conts,col=adjustcolor("black",0.25))
    points(tort.coords, pch=20, cex=1)

# look at with tortoise whatnot
layout(matrix(1:6,nrow=2))
for (k in 1:5) {
    plot_wbd(k, col=cols[[k]], labels=NULL)
    lines(conts,col=adjustcolor("black",0.25))
    points(tort.coords, pch=20, cex=3)
}
