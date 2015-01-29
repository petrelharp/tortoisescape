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

## and the north-south masks

nus <- raster("nussear.grd")
masked <- raster("mask_crew_dem_2K_sea.grd")

lon <- raster("lon_utm_30.grd")
lat <- raster("lat_utm_30.grd")

# north mask
n.mask <- mask( masked, (lon>5e5) & (lat>3.8e6), maskvalue=FALSE )
n.mask.clump <- clump(n.mask)
big.clump <- which.max(table(values(n.mask.clump)))
n.mask[n.mask.clump != big.clump ] <- NA
writeRaster(n.mask,file="mask_crew_dem_2K_sea_north.grd",overwrite=TRUE)

# south mask
s.mask <- mask( masked, (lat>4e6) & (lon>4.2e5), maskvalue=TRUE )
s.mask.clump <- clump(s.mask)
big.clump <- which.max(table(values(s.mask.clump)))
s.mask[s.mask.clump != big.clump ] <- NA
writeRaster(s.mask,file="mask_crew_dem_2K_sea_south.grd",overwrite=TRUE)
