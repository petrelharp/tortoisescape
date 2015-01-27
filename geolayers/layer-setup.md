We have a new set of layers, and need to set up inference on them.
These are in `geolayers/expanded/expanded-TIFF/` (not really, but symlinked in).
First, we need to make the scaled down layers: in the directory `geolayers/`, run:
```
./aggregate-batch.sh expanded/expanded-TIFF/ expanded/
```
This calls the script `aggregate-layers.R` on each layer in the directory `expanded/expanded-TIFF/`,
aggregating by a factor of 2, and repeating.

Then, we need to put together the consensus masking NA layer (this is in `make-masking-layer.R`):
```
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
```

Hm, still gotta do some processing to get the consensus NA layer.
```
require(raster)
crew <- raster("crew_30_NA")
dem.mask <- raster("dem_2K_mask.tif")
water <- raster("water0_30.tif")
sea <- raster("sea_NA_30.tif")
masked <- ( is.na(crew) | is.na(dem.mask) | !is.na(sea) )
masked[!masked] <- NA
```

And, let's smooth out the edges some, extending it out by 300m:
```
d.mask <- distance(masked)
mask.exp <- ( d.mask > 300 )
mask.exp[mask.exp==0] <- NA
```
and by removing disconnected pieces:
```
mask.clumps <- clump(mask.exp)
clump.sizes <- table(values(mask.clumps))
big.clump <- which.max( clump.sizes )
mask.exp[mask.clumps!=big.clump] <- NA
```



