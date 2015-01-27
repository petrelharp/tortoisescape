#!/usr/bin/Rscript

usage <- "Combine crew, dem_2K_mask, and sea_NA_30 into a single layer for masking.
Usage:
    Rscript make-masking-layer.R (name of directory) (aggregation factor for clumping)
"

argvec <- if (interactive()) { scan(what='char') } else { commandArgs(TRUE) }
if ( length(argvec)<2 ) { stop(usage) }

require(raster)

subdir <- argvec[1]
ag.fact <- as.numeric(argvec[2])

setwd(subdir)
grd.or.tif <- function (x) { paste(x, if (file.exists(paste(x,".grd",sep=''))) { ".grd" } else { ".tif" }, sep='' ) }
crew <- raster(grd.or.tif("crew_30_NA"))
dem.mask <- raster(grd.or.tif("dem_2K_mask"))
sea <- raster(grd.or.tif("sea_NA_30"))
masked <- mask(crew,dem.mask)
masked <- mask(masked,sea)
inv.masked <- masked
inv.masked[is.na(masked)] <- 1
inv.masked[!is.na(masked)] <- NA
d.mask <- distance(inv.masked)
mask.exp <- ( d.mask > 300 )
mask.exp[mask.exp==0] <- NA
mask.ag <- if (ag.fact>1) { aggregate(mask.exp,fact=ag.fact,fun=mean,na.rm=FALSE) } else { mask.exp }
mask.clumps <- clump(mask.ag)
clump.sizes <- table(values(mask.clumps))
big.clump <- which.max( clump.sizes )
mask.ag[mask.clumps!=big.clump] <- NA
mask.disag <- if (ag.fact>1) { disaggregate( mask.ag, fact=ag.fact ) } else { mask.ag }
masked <- mask( masked, crop(mask.disag,masked) )
writeRaster(masked,file="mask_crew_dem_2K_sea.grd")
removeTmpFiles()
