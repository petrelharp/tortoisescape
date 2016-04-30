library(raster)

load("county_lines.Robj")  # provides county_lines
dem <- raster("dem_30")    # raster of elevation

# the PROJ.4 string used by our rasters
# note this is:
#     "+proj=eqdc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
.raster.proj4 <- CRS(scan("../raster_GCS_CRS_proj4.txt",what="char",sep="\n"))
.raster.crs <- CRS(our.proj4)
