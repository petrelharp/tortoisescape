require(raster)
setwd("TIFFs")
crew <- raster("crew_NA_1km.tif")
dem <- raster("dem_2k_1km.tif")
sea <- raster("sea_NA_1km.tif")
nus <- raster("nussear.tif")
crew <- crop(crew,nus)
dem <- crop(dem,nus)
sea <- crop(sea,nus)
nus <- crop(nus,sea)
masked <- mask(crew,dem)
masked <- mask(masked,sea)
masked <- mask(masked,nus)
# remove isolated regions
mask.clumps <- clump(masked,directions=4)
clump.sizes <- table(values(mask.clumps))
big.clump <- which.max( clump.sizes )
masked[mask.clumps!=big.clump] <- NA
writeRaster(masked,file="mask_crew_dem_2K_sea.tif",overwrite=TRUE)

setwd("..")
writeRaster(crew,"crew_NA_1km.grd",overwrite=TRUE)
writeRaster(dem,"dem_NA_1km.grd",overwrite=TRUE)
writeRaster(sea,"sea_NA_1km.grd",overwrite=TRUE)
writeRaster(nus,"nussear.grd",overwrite=TRUE)
writeRaster(masked,file="mask_crew_dem_2K_sea.grd",overwrite=TRUE)

removeTmpFiles()


