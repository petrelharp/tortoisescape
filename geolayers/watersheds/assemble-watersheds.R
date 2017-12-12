
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
sessionInfo()
# counties <- get_counties(dem)
# conts <- get_contours(dem)
# shade <- get_shading(dem)
print(getwd())
num.wbd <- 3 # note this is not doing WBD 12, which needs more RAM
wbd.list <- lapply( 2*(1:num.wbd), function (wbd_num) {
            gClip( 
                  spTransform( 
                               readOGR("WBD/WBD.gdb",paste0("WBDHU",wbd_num)),
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
#     plot(shade, col=adjustcolor(grey(seq(0,1,length.out=101)),0.5), legend=FALSE, add=TRUE )
#     lines(conts,col=adjustcolor("black",0.25))
#     lines(counties, lwd=0.5, col=adjustcolor("red",0.5))
    if (do.torts) points(tort.coords, pch=20, cex=1)
    xy <- gCentroid( wbd.list[[k]], byid=TRUE )
    text(xy,labels)
}

pdf(file="watersheds.pdf", width=8, height=8, pointsize=10)
for (k in 1:num.wbd) {
    plot_wbd(k, labels=NULL, main=paste("WBD",2*k))
}
dev.off()

# NB not adapted for num.wbd variable!
if (FALSE) {
    # check for contiguity
    layout(matrix(1:8,nrow=2))
    for (k in 2:5) {
        plot_wbd(k-1, do.torts=FALSE)
        plot_wbd(k, labels=hierarchies[[k]], do.torts=FALSE)
    }


    layout(t(1:2))
    plot_wbd(4, labels=NULL )
    plot_wbd(5, labels=NULL )
}
