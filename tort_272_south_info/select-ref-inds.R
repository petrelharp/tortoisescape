#!/usr/bin/Rscript

require(raster)
require(rgdal)

source("../inference/resistance-fns.R")

layer <- raster("../visualization/crm_dem_30.gri")
sample.locs.obj <- load("geog_coords.RData")
sample.locs <- spTransform(get(sample.locs.obj),CRSobj=CRS(proj4string(layer)))

# 8.5km chosen to get 47 of these
nonmissing <- which(!is.na(values(layer)))
neighborhoods <- get.neighborhoods(8500, sample.locs, nonmissing, layer, numcores=getcores(), na.rm=TRUE )
nonoverlapping <- which.nonoverlapping(neighborhoods)

ref.names <- row.names(sample.locs)[nonoverlapping]

cat(paste('"',paste(ref.names,collapse='", "'),'"\n',sep=''))
