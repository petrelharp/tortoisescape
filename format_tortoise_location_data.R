################################################################
################################################################
#	get raster data for all tortoise locations
################################################################
################################################################

if(file.exists("~/Desktop/Dropbox/tortoisescape")){
	setwd("~/Desktop/Dropbox/tortoisescape")
}

require(raster)
require(rgdal)

rasterOptions(tmpdir=".")

################################
#	getting tortoise locations
################################

#	First, read in tortoise data from .csv file
tort.coord.data <- read.csv("1st_180_torts.csv")

#	also, grab the proj4 string describing the GCS of
#		the rasters from jannet
load("raster_GCS_CRS_proj4.Robj")
tort.coords <- cbind(tort.coord.data$Easting,tort.coord.data$Northing)

#	get the UTM zone, because not all are from zone 11
tort.zone <- tort.coord.data$UTM_Zone

#	go through tortoise coordinates and,
#		using the appropriate UTM zone,
#		convert them into SpatialPoints objects
#		and then change their proj4 string to that
#		of the rasters.
tort.coords_zone11 <- SpatialPoints(tort.coords[grepl("11",tort.zone),],proj4string=CRS("+proj=utm +zone=11"))
tort.coords_zone12 <- SpatialPoints(tort.coords[grepl("12",tort.zone),,drop=FALSE],proj4string=CRS("+proj=utm +zone=12"))
tort.coords_zone11_rasterGCS <- spTransform(tort.coords_zone11,raster_GCS_CRS_proj4)
tort.coords_zone12_rasterGCS <- spTransform(tort.coords_zone12,raster_GCS_CRS_proj4)
tort.coords.rasterGCS <- rbind(tort.coords_zone11_rasterGCS[1:(grep("12",tort.zone)-1),],
									tort.coords_zone12_rasterGCS,
								tort.coords_zone11_rasterGCS[grep("12",tort.zone):(length(tort.zone)-1),])

row.names(tort.coords.rasterGCS) <- tort.coord.data$EM_Tort_ID
save(tort.coords.rasterGCS,file="tort.coords.rasterGCS.Robj")

################################
#	making pairwise distance table
################################

tort.distance <- pointDistance(tort.coords.rasterGCS,lonlat=FALSE)


