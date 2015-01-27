#!/usr/bin/Rscript

usage <- "Combine crew, dem_2K_mask, and sea_NA_30 into a single layer for masking.
Usage:
    Rscript make-masking-layer.R (name of directory) (aggregation factor for clumping)
"

argvec <- if (interactive()) { scan(what='char') } else { commandArgs(TRUE) }
if ( length(argvec)<2 ) { stop(usage) }

require(raster)
rasterOptions(tmpdir=".")

subdir <- argvec[1]
ag.fact <- argvec[2]

setwd(subdir)
grd.or.tif <- function (x) { paste(x, if (file.exists(paste(x,".grd",sep=''))) { ".grd" } else { ".tif" }, sep='' ) }
crew <- raster(grd.or.tif("crew_30_NA"))
dem.mask <- raster(grd.or.tif("dem_2K_mask"))
water <- raster(grd.or.tif("water0_30"))
sea <- raster(grd.or.tif("sea_NA_30"))
masked <- mask(crew,dem.mask)
masked <- mask(masked,sea)
d.mask <- distance(masked)
mask.exp <- ( d.mask > 300 )
mask.exp[mask.exp==0] <- NA
mask.ag <- aggregate(mask.exp,fact=ag.fact,fun=mean,na.rm=FALSE)
mask.clumps <- clump(mask.ag)
clump.sizes <- table(values(mask.clumps))
big.clump <- which.max( clump.sizes )
mask.ag[mask.clumps!=big.clump] <- NA
mask.disag <- disaggregate( mask.ag, fact=ag.fact )
masked <- mask( masked, mask.disag )
writeRaster(masked,file="mask_crew_dem_2K_sea.grd")
removeTmpFiles()
