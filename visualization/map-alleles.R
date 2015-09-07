#!env Rscript

usage <- "
Make maps for each allele in files produced by extract-random-alelles.R:
    Rscript map-random-alleles.R [list of directories to do this in]
"

dirnames <- if (interactive()) { scan(what='char') } else { commandArgs(TRUE) }

if (length(dirnames)==0) { stop(usage) }

####### make pngs
# set up for mapping
coord.obj <- load(file.path(tortdir,"tort_272_info","geog_coords.RData"))
coords <- get(coord.obj)
tort.ids <- row.names(coords)
bases <- c("A","C","G","T")
base.cols <- adjustcolor(1:4,.5)
nsamples <- length(count.header)/4
stopifnot( length(tort.ids) == nsamples )

library(raster)
library(maptools)
layer <- raster(file.path(tortdir,"visualization/dem_30"))
player <- function (main='') { plot(layer,legend=FALSE,xlab="",ylab="",xaxt="n",yaxt="n",legend.mar=0,box=FALSE,main=main) }
# read in other info
pcs <- read.csv(file.path(tortdir,"tort_272_info","pcs.csv"),header=TRUE,stringsAsFactors=FALSE)
stopifnot( all( tort.ids %in% pcs$etort ) )
pc.cols <- adjustcolor( ifelse( pcs$PC1[match(tort.ids,pcs$etort)] > 0, "blue", "purple" ), .75 )



for (workdir in dirnames) {
    site.info <- read.table( file.path(workdir, "random_sites.info"), header=TRUE )
    counts <- as.matrix( read.table( file.path(workdir, "random_sites.counts"), header=TRUE ) )
    for (site in 1:nrow(site.info)) {
        png(
            file=file.path(workdir,paste(site.info$chr[site],"_",site.info$pos[site],".png",sep="")), 
            width=4*144, height=4*144, pointsize=10, res=144 )
        player(main=paste(site.info$chr[site],site.info$pos[site],"depth:"site.info$totDepth[site],"freq:",site.info$knownEM[site],"nInd:",site.info$nInd[site]))
        for (k in seq_along(bases)) {
            nn <- k+4*(0:(nsamples-1))
            points(coords,pch=20,cex=counts[site,nn],col=base.cols[k])
        }
        dev.off()
    }
}
