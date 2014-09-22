# script to build an adjacency matrix accounting for the missing values

require(Matrix)
require(raster)

# load raster (change to full size for actual use)
rast<- raster("../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_m2_02tmax")
## For cluster
# setwd("/home/rcf-40/pralph/panfs/tortoises/tortoisescape")
# layer.prefix <- c("geolayers/TIFF/masked/")
# layer.names <- list.files( layer.prefix, "*.grd" )
# rast <- raster(paste(layer.prefix,layer.names[1],sep=''))

layer <- values(rast)

# find dimensions of grid
n <- dim(rast)[2]; m <- dim(rast)[1]

# find locations of entriea
loc <- which(!is.na(layer))

# clear "layer" from memory
#rm(layer)
#gc()

# use function "adjacent" (returns )
system.time(ij <- adjacent(rast,cells=loc,target=loc)) # to and from cells both loc

# clear rast and loc from memory
rm(rast,loc)
gc()

# declare matrix
system.time(G <- sparseMatrix(ij[,1],ij[,2],x=rnorm(length(ij[,1]),0,1)))
