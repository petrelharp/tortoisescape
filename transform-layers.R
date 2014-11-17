#!/usr/bin/Rscript

# make some new layers
# that are transforms of the old

require(raster)
rasterOptions(tmpdir=".")

basedir <- "geolayers/TIFF/"

layer.name <- "dem_30"
new.name <- "dem_30_m800_sq"
xform <- function (x) { scale( (x-800)^2 ) }
maxval <- 3

rast <- raster( paste(basedir,"masked/crop_resampled_masked_",layer.name,sep='') )
new.rast <- xform(rast)
new.rast[new.rast>maxval] <- maxval
writeRaster(new.rast, file=paste(basedir,"masked/crop_resampled_masked_",new.name,sep=''), overwrite=TRUE )


# from crop_mask_rasters.R
for (fact in c(10,100,500)) {
    prefix <- if (fact %in% c(10,100)) { "crop_resampled_masked_aggregated_" } else { "" }
    outfile <- paste(basedir,fact,"x/",prefix,fact,"x_",new.name,sep='')
    aggregate(new.rast,
                fact=fact,
                fun=mean,
                na.rm=TRUE,
                filename=outfile,
                overwrite=TRUE)
}

