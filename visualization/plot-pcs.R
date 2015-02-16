#!/usr/bin/Rscript

usage <- "Makes various plots of PCs.
Usage:
    Rscript plot-pcs.R (.csv with pcs)
which will then make several image files in the same directory the file with pcs is.
"

argvec <- if (interactive()) { scan(what="char") } else { commandArgs(TRUE) }

if (length(argvec)<1) { stop(usage) }
infile <- argvec[1]
outbase <- file.path( dirname(infile), gsub("[.].*$","-",basename(infile)) )
pcs <- read.csv(infile)

require(raster)
require(colorspace)
require(maps)
require(maptools)
dem <- raster("../visualization/dem_30.gri")
load("../visualization/county_lines.Robj")  # provides county_lines
sample.loc.obj <- load("geog_coords.RData")
sample.locs <- get(sample.loc.obj)
sample.coords <- coordinates(sample.locs)

colorize <- function (x,alpha=1) { diverge_hcl(128, h=c(225,0), c=100, l=c(60,95), alpha=alpha)[cut(x,128)] }

plot.vec <- function (x,...) {
    plot(dem,legend=FALSE,...)
    lines(county_lines)
    x.fac <- cut(x,128)
    cols <- diverge_hcl(128, h=c(225,0), c=100, l=c(60,95))
    points(sample.locs,pch=21,cex=2,col=grey(.2), bg=cols[as.numeric(x.fac)] )
    tmp <- dem; tmp[] <- NA
    tmp[1:length(unique(cols))] <- tapply( x, x.fac, mean )
    plot( tmp, legend.only=TRUE, legend.width=2,
            breaks=pretty(x,8),
            col=cols[(1:8)*length(cols)/8]
        )
}


# plots of PCs against each other, colored by northing and easting
outfile <- paste(outbase,"PC-maps.pdf",sep='')
cat("Plotting to ", outfile, "\n")
pdf(file=outfile, width=10, height=8, pointsize=10)
    pairs( pcs[c("PC1","PC2","PC3","PC4")], 
        lower.panel=function(x,y,...){points(x,y,bg=colorize(sample.coords[,2],.5), pch=21, cex=2, col=grey(.20) )}, 
        upper.panel=function(x,y,...){points(x,y,bg=colorize(sample.coords[,1],.5), pch=21, cex=2, col=grey(.20) )},
        main="Northing (above), Easting (below)" )
dev.off()

# maps, with PCs on
outfile <- paste(outbase,"maps-with-PCs.pdf",sep='')
cat("Plotting to ", outfile, "\n")
pdf(file=outfile, width=10, height=8, pointsize=10)
    plot(dem,main="Elevation with tortoise IDs")
    lines(county_lines)
    text(sample.locs,labels=gsub("etort.","",row.names(sample.locs)))
    ncols <- 16
    for (pc in c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")) {
        plot.vec(pcs[[pc]],main=pc)
    }
dev.off()


