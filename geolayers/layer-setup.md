We have a new set of layers, and need to set up inference on them.
These are in `geolayers/expanded/expanded-TIFF/` (not really, but symlinked in).
First, we need to make the scaled down layers: in the directory `geolayers/`, run:
```
./aggregate-batch.sh expanded/expanded-TIFF/ expanded/
```
This calls the script `aggregate-layers.R` on each layer in the directory `expanded/expanded-TIFF/`,
aggregating by a factor of 2, and repeating.

Hm, still gotta do some processing to get the consensus NA layer.
```
require(raster)
rasterOptions(tmpdir=".")
crew <- raster("crew_30.tif")
dem.mask <- raster("dem_2K_mask.tif")
water <- raster("water0_30.tif")
sea <- raster("sea_NA_30.tif")
combined.mask <- ( (crew>0) & !is.na(dem.mask) & !is.na(sea) )
```

And, let's smooth out the edges some, extending it out by 300m:
```
masked <- combined.mask
masked[!combined.mask] <- NA
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



