require(raster)
setwd("TIFFs")
crew <- raster("crew_NA_1km.tif")
dem.mask <- raster("dem_2k_1km.tif")
sea <- raster("sea_NA_1km.tif")
masked <- mask(crew,dem.mask)
masked <- mask(masked,sea)
mask.clumps <- clump(masked)
clump.sizes <- table(values(mask.clumps))
big.clump <- which.max( clump.sizes )
mask.ag[mask.clumps!=big.clump] <- NA
masked <- mask( masked, mask.disag )
writeRaster(masked,file="mask_crew_dem_2K_sea.tif")

setwd("..")
tiffs <- list.files("TIFFs","tif$",full.names=TRUE)
grds <- gsub("tif$","grd",basename(tiffs))
for (k in seq_along(tiffs)) { r <- raster(tiffs[k]); writeRaster(r,file=grds[k]) }

removeTmpFiles()


