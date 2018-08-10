# Modified from locate_tortoises.R
#   Run this in e.g. tort_272_info

library(rgdal)
library(raster)
library(sp)
library(rgeos)
library(spdep)
library(colorspace)
library(maps)
source("../visualization/map-utils.R",chdir=TRUE)
watershed_dir <- "../geolayers/watersheds"
source(file.path(watershed_dir, "watersheds-fns.R"), chdir=TRUE)

tort.coords.obj <- load("geog_coords.RData")
tort.coords <- get(tort.coords.obj)

dem <- get_elev()
conts <- get_contours(dem)
shade <- get_shading(dem)
counties <- get_counties(dem)
ocean <- get_ocean(dem)

load("handpicked_WBD8.RData")
tort_wbd <- read.csv("watershed_assignments.csv", header=TRUE)

# change this to change the region colors
# region_colors <- terrain.colors(length(sub_wbd8))
region_colors <- adjustcolor(rainbow(length(sub_wbd8)), 0.75)
# and the labels
region_labels <- seq_along(sub_wbd8)

# png(file="fancy_watershed_assignments.png", width=2.75*288, height=3*288, pointsize=10, res=288)
pdf(file="fancy_watershed_assignments.pdf", width=3, height=3, pointsize=10)
par(mar=c(1.0,1.0,1.0,1.0))
    plot(sub_wbd8, col=adjustcolor(region_colors, 0.5), lwd=1, 
             xaxt='n', yaxt='n', xlab='', ylab='',
             xlim=.expand(extent(tort.coords)[1:2],.15), 
             ylim=.expand(extent(tort.coords)[3:4],.15))
    plot(shade, col=adjustcolor(grey(seq(0,1,length.out=101)),1.0), legend=FALSE, 
         add=TRUE )
    plot(ocean, add = TRUE, col = "light blue")
    # lines(conts,col=adjustcolor("black",0.25))
    lines(counties, lwd=0.5, col=adjustcolor("red",0.5))
    points(tort.coords, pch=20, cex=1)
    xy <- gCentroid( sub_wbd8, byid=TRUE )
    points(xy, cex=3, pch=20, col=adjustcolor(region_colors, 0.5))
    text(xy, labels, cex=1, col='black')
    scalebar(5e4, label=c(NA,"50km",NA))
dev.off()
