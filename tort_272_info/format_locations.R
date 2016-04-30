
library(sp)
library(rgdal)

torts <- read.csv("272torts_metadata.csv",stringsAsFactors=FALSE)

#	this is the proj4 string describing the GCS of the rasters from jannet
raster.proj4 <- CRS(scan("../raster_GCS_CRS_proj4.txt",what="char",sep="\n"))

# obtained PROJ.4 string from http://spatialreference.org/ref/epsg/2153/proj4/
nad.11 <- spTransform( SpatialPoints(torts[,c("Easting","Northing")], proj4string=CRS("+proj=utm +zone=11 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")), raster.proj4)
# obtained PROJ.4 string from http://spatialreference.org/ref/epsg/wgs-84-utm-zone-11n/
wgs.11 <- spTransform( SpatialPoints(torts[,c("Easting","Northing")],proj4string=CRS("+proj=utm +zone=11 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")), raster.proj4)

stopifnot( all( abs(pointDistance(nad.11,wgs.11)) < 1e-3 ) )

tort.coords <- nad.11
row.names(tort.coords) <- torts$EM_Tort_ID
save(tort.coords,file="geog_coords.RData")

# get long/lat:

tort.coords.longlat <- cbind( coordinates( spTransform( tort.coords, CRS("+proj=longlat") ) ), torts$EM_Tort_ID )
colnames(tort.coords.longlat) <- c("longitude","latitude","id")
write.csv(tort.coords.longlat, row.names=FALSE, file="long-lat.csv")

# get geographic distances

dists <- pointDistance( tort.coords, lonlat=FALSE,  )
ut <- upper.tri(dists,diag=TRUE)
dists.df <- data.frame(
                       etort1=row.names(tort.coords)[row(dists)][ut],
                       etort2=row.names(tort.coords)[col(dists)][ut],
                       distance=dists[ut] )

write.csv(dists.df, row.names=FALSE, file="geog_distance.csv")

if (FALSE) {
    ## data checking, correcting, and exploring

    # these look reasonable
    layout(t(1:2))
    plot(dem, main="zone 11")
    lines(county_lines)
    points(nad.11, pch=1)
    points(wgs.11, pch=1, col='red')

    # Clearly, all are in zone 11.
    if (FALSE) {
        nad.12 <- spTransform( SpatialPoints(torts[,c("Easting","Northing")],proj4string=CRS("+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")), raster.proj4)
        wgs.12 <- spTransform( SpatialPoints(torts[,c("Easting","Northing")],proj4string=CRS("+proj=utm +zone=12 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")), raster.proj4)
        plot(dem, main="zone 12")
        lines(county_lines)
        points(nad.12, pch=2)
        points(wgs.12, pch=2, col='red')
    }

    # how far apart are they?
    nad.wgs.dist <- pointDistance(nad.11,wgs.11)
    # all within 0.0001098557
    range(nad.wgs.dist)

    # fix up etort-231: this is now FIXED in 272torts_metadata.csv
    # these are next to etort-231:
    points(nad.11[c(216,240,262)],col='green',pch=20)
    new.231 <- spTransform( SpatialPoints( cbind( Easting=711975, Northing=3823493), proj4string=CRS("+proj=utm +zone=11 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")), raster.proj4)
    points(new.231, col='blue', cex=2, lwd=2)

}
