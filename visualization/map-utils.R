library(raster)
library(maps)
library(maptools)

# the PROJ.4 string used by our rasters
# note this is:
#     "+proj=eqdc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
.raster.proj4 <- scan("../raster_GCS_CRS_proj4.txt",what="char",sep="\n")
.raster.crs <- CRS(.raster.proj4)

.thisdir <- file.path(normalizePath("."))

.ivanpah.bbox <- structure(c(-1754578, -233915.5, -1722266, -176116.9),
    .Dim = c(2L, 2L), .Dimnames = list(c("Easting", 
    "Northing"), c("min", "max")))

#' Get county outlines to overlay on a spatial object
#' @param x The spatial object
get_counties <- function (x) {
    the.proj4 <- if (missing(x)) { .raster.proj4 } else { proj4string(x) }
    mojave.outlines <- map("county",c("california","nevada","arizona"),
        xlim=c(-120,-113),ylim=c(31,37), plot=FALSE)
    mojave.outlines.sp <- map2SpatialLines(mojave.outlines,proj4string=CRS("+proj=longlat"))
    spTransform(mojave.outlines.sp,CRS(the.proj4))
}

#' Get state lines to overlay on a spatial object
#' @param x The spatial object
get_statelines <- function (x) {
    the.proj4 <- if (missing(x)) { .raster.proj4 } else { proj4string(x) }
    state.lines <- map(database="state",regions=c("California","Arizona","Nevada"),plot=FALSE)
    state.lines.spobj <- map2SpatialLines(state.lines,proj4string=CRS("+proj=longlat"))
    spTransform(state.lines.spobj,CRS(the.proj4))
}


#' Get an elevation raster matching with a spatial object
#' @param x The spatial object
get_elev <- function (x) {
    dem <- raster(file.path(.thisdir,"dem_30"))
    the.proj4 <- if (missing(x)) { .raster.proj4 } else { proj4string(x) }
    if (the.proj4 != proj4string(dem)) {
        dem <- spTransform( rasterToContour(dem,nlevels=25), CRS(the.proj4) )
    }
    return(dem)
}

#' Get elevation contours to overlay on a spatial object
#' @param x The spatial object
get_contours <- function (x) {
    conts <- rasterToContour( raster(file.path(.thisdir,"dem_30")), nlevels=25 )
    the.proj4 <- if (missing(x)) { .raster.proj4 } else { proj4string(x) }
    if (the.proj4 != proj4string(conts)) {
        conts <- spTransform( conts, CRS(the.proj4) )
    }
    return( conts )
}

#' Get elevation shading to overlay on a spatial object
#' @param x The spatial object
#' Plot the result with, e.g.
#' plot( shade, col=adjustcolor(grey(seq(0,1,length.out=101)),0.5), legend=FALSE )
get_shading <- function (x) {
    SR <- raster(file.path(.thisdir,"cropped_SR_HR.tif"))
    the.proj4 <- proj4string(x)
    projectRaster(SR,to=raster(extent(x),res=500,crs=CRS(the.proj4)))
}

#' Get ocean to overlay
#' @param x The spatial object
#' plot(ocean, add = TRUE, col = "light blue")
get_ocean <- function (x) {
    the.proj4 <- if (missing(x)) { .raster.proj4 } else { proj4string(x) }
    xbox <- as(extent(x),"SpatialPolygons")
    crs(xbox) <- CRS( the.proj4 )
    if (FALSE) {
        # this is hella messy; gotta clean it up
        land <- spTransform( readOGR(file.path(.thisdir,"natural_earth"),"ne_10m_land"), crs(the.proj4) )
        tmp <- land
        # this takes FOREVER
        tmp@polygons <- lapply(land@polygons, checkPolygonsHoles)
        tmp2 <- gBuffer(tmp, byid=TRUE, width=0)
        bigbox <- SpatialPolygons(
                    list(Polygons(list(
                           Polygon(cbind(c(-4e6,-1e6)[c(1,1,2,2)],
                                         c(-8e5,4e5)[c(1,2,2,1)]))), ID='box')),
                    proj4string=CRS(the.proj4))
        land <- crop(tmp2, bigbox)
        save(land, file="natural_earth_cleaned.RData")
    }
    load(file.path(.thisdir, "natural_earth_cleaned.RData"))
    landx <- rgeos::gDifference(xbox, crop(land, xbox))
    return(landx)
}

