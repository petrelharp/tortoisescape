#!/usr/bin/Rscript
require(raster)
require(rgdal)

pref <- raster("pref_alt_poly.tif")
alt1 <- raster("alt_1_poly.tif")
alt2 <- raster("alt_2_poly.tif")
alt3 <- raster("alt_3_poly.tif")
alt4 <- raster("alt_4_poly.tif")

