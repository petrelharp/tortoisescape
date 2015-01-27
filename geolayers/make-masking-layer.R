#!/usr/bin/Rscript

usage <- "Combine crew, dem_2K_mask, and sea_NA_30 into a single layer for masking.
Usage:
    Rscript make-masking-layer.R (name of directory) ( more directories )
"

if ( length(commandArgs(TRUE))<1 ) { stop(usage) }

require(raster)
rasterOptions(tmpdir=".")

orig.dir <- getwd()

subdir <- commandArgs(TRUE)[1]
setwd(subdir)
crew <- raster("crew_30_NA")
dem.mask <- raster("dem_2K_mask")
water <- raster("water0_30")
sea <- raster("sea_NA_30")
masked <- mask(crew,dem.mask)
masked <- mask(masked,sea)
d.mask <- distance(masked)
mask.exp <- ( d.mask > 300 )
mask.exp[mask.exp==0] <- NA
mask.ag <- aggregate(mask.exp,fact=64,fun=mean,na.rm=FALSE)
mask.clumps <- clump(mask.ag)
clump.sizes <- table(values(mask.clumps))
big.clump <- which.max( clump.sizes )
mask.ag[mask.clumps!=big.clump] <- NA
mask.disag <- disaggregate( mask.ag, fact=64 )
masked <- mask( masked, mask.disag )
writeRaster(masked,file="mask_crew_dem_2K_sea.grd")
