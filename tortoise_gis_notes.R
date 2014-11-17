################################################################
################################################################
#	NOTES ON WORKING WITH GIS LAYERS IN R	
################################################################
################################################################
require(raster)
# this put the temp stuff in my temp directory, and filled up that partition:
rasterOptions(tmpdir=".")

# on gideon's computer, which is too small to 
#	have the raster tifs on HD
if(file.exists("/Volumes/BBURD/tortTestGIS")){
	setwd("/Volumes/BBURD/tortTestGIS")
}

if(file.exists("/Volumes/cooplab1/tortoises/geolayers/10x")){
	setwd("/Volumes/cooplab1/tortoises/geolayers/10x")
}

# go through all rasters and get the value range for each one
raster.list <- list.files()[grep("*.gri",list.files(pattern="_10x"))]
raster.list.values <- vector("list",length=length(raster.list))
for(i in 1:length(raster.list)){
	tmp.raster <- raster(raster.list[i])
	raster.list.values[[i]] <- list(min = tmp.raster@data@min, 
									max = tmp.raster@data@max,
									nodatavalue = tmp.raster@file@nodatavalue)
}
raster.names <- gsub("crop_resampled_masked_aggregated_10x_","",raster.list)
	raster.names <- gsub(".gri","",raster.names)
	raster.names <- gsub("crop_masked_aggregated_10x_","",raster.names)
names(raster.list.values) <- raster.names
raster.value.dataframe <- data.frame("min.value" = unlist(lapply(raster.list.values,"[[",1)),
									"max.value" = unlist(lapply(raster.list.values,"[[",2)),
									"no.data.value" = unlist(lapply(raster.list.values,"[[",3)),
									row.names=raster.names)
save(raster.value.dataframe,file="raster.value.dataframe.Robj")

	
	
# everything below in the if(FALSE) statement
#	was exploratory code I made for working
#	with rasters.  Doesn't have to be run
#	for anything now, but good to keep old notes around

if(FALSE){
# raster is clever about what it chooses to read into working
#	memory.  this is great for big data files, especially
#	compared to the raw readGDAL() function.  raster also has
#	nice vignettes: 
#		vignette("Raster",package="raster")
#		vignette("functions",package="raster")
#		vignette("rasterfile",package="raster")

# create a raster of one of the .tifs.  this creates a link
#	between R and the specific file, from which it will dynamically
#	pull info required of it.
road <- raster("geolayers/TIFF/road_30.tif")
precip <- raster("geolayers/TIFF/annual_precip.tif")
aspect <- raster("geolayers/TIFF/aspect_30.tif")
lat <- raster("geolayers/TIFF/lat_gcs_30.tif")
long <- raster("geolayers/TIFF/lon_gcs_30.tif")

extent(road)
extent(precip)
extent(aspect)
extent(lat)
extent(long)

# lat and long and precip have the same extent
# road and aspect have the same extent
# So, in the next lines of code, we want to get all these 5 rasters to the same extent
road@extent@xmin	#-2154584
lat@extent@xmin		#-2154993	#lat goes farther west

road@extent@xmax	#-1497096	#lat goes farther east
lat@extent@xmax		#-1475779

road@extent@ymin	#-627324.3	#they both go equally south
lat@extent@ymin		#-627324.3

road@extent@ymax	#33770.03	#lat goes farther north
lat@extent@ymax		#41742.17


#	so, we first crop lat, long, and precip to road
writeRaster(
			crop(lat,
				extent(
					max(road@extent@xmin,
							lat@extent@xmin),
					min(road@extent@xmax,
							lat@extent@xmax),
					max(road@extent@ymin,
							lat@extent@ymin),
					min(road@extent@ymax,
							lat@extent@ymax))),
				filename="crop_lat",overwrite=TRUE)

writeRaster(
			crop(long,
				extent(
					max(road@extent@xmin,
							lat@extent@xmin),
					min(road@extent@xmax,
							lat@extent@xmax),
					max(road@extent@ymin,
							lat@extent@ymin),
					min(road@extent@ymax,
							lat@extent@ymax))),
				filename="crop_long",overwrite=TRUE)				

writeRaster(
			crop(precip,
				extent(
					max(road@extent@xmin,
							lat@extent@xmin),
					min(road@extent@xmax,
							lat@extent@xmax),
					max(road@extent@ymin,
							lat@extent@ymin),
					min(road@extent@ymax,
							lat@extent@ymax))),
				filename="crop_precip",overwrite=TRUE)

cropped_lat <- raster("crop_lat.grd")
cropped_long <- raster("crop_long.grd")
cropped_aspect <- raster("crop_precip.grd")

#	then we crop road to the cropped lat (because, for whatever reason, crop() overshoots it a bit
writeRaster(
			crop(road,
				extent(
					max(road@extent@xmin,
							cropped_lat@extent@xmin),
					min(road@extent@xmax,
							cropped_lat@extent@xmax),
					max(road@extent@ymin,
							cropped_lat@extent@ymin),
					min(road@extent@ymax,
							cropped_lat@extent@ymax))),
				filename="crop_road",overwrite=TRUE)

writeRaster(
			crop(aspect,
				extent(
					max(road@extent@xmin,
							cropped_lat@extent@xmin),
					min(road@extent@xmax,
							cropped_lat@extent@xmax),
					max(road@extent@ymin,
							cropped_lat@extent@ymin),
					min(road@extent@ymax,
							cropped_lat@extent@ymax))),
				filename="crop_aspect",overwrite=TRUE)

cropped_aspect <- raster("crop_aspect.grd")
cropped_road <- raster("crop_road.grd")

#	Then we resample road and aspect to lat
resample(cropped_road,cropped_lat,filename="resampled_cropped_road",overwrite=TRUE)
resample(cropped_aspect,cropped_lat,filename="resampled_cropped_aspect",overwrite=TRUE)

resampled_cropped_road <- raster("resampled_cropped_road")
resampled_cropped_aspect <- raster("resampled_cropped_aspect")

#	and now:
compareRaster(resampled_cropped_road,cropped_lat)
compareRaster(resampled_cropped_aspect,cropped_lat)
}


# the code below explores how
#	to get the tortoise coordinates
#	into the coordinate space of rasters
#	we have from jannet
require(rgdal)
require(maps)
require(maptools)
require(sp)

#	read in some sample rasters
resampled_cropped_road <- raster("resampled_cropped_road")
resampled_cropped_aspect <- raster("resampled_cropped_aspect")
cropped_lat <- raster("crop_lat.grd")
cropped_long <- raster("crop_long.grd")
cropped_precip <- raster("crop_precip.grd")

#	read in the tortoise metadata for the first
#		180 sequenced tortoises
tort.coords <- read.csv("1st_180_torts.csv")

#	grab their lat/long, which are in UTM
utm.coords <- cbind(tort.coords$Easting,tort.coords$Northing)

#	make those into an object of class SpatialPoints,
#		using the proj4string for utm coords in zone 11 (where they are)
utm.coords <- SpatialPoints(utm.coords,proj4string=CRS("+proj=utm +zone=11"))

#	project those UTM coordinates into lat and long
#		(just for easy plotting, and as a sanity check)
latlong.coords <- spTransform(utm.coords,CRS("+proj=longlat +init:epsg:2955"))

#	project those UTM coordinates into the reference
#		system of the rasters from Jannet.
#	note that it seems we can do this just pulling out the proj4string
#		from the raster object.
proj.latlong.coords <- spTransform(utm.coords,cropped_lat@crs)


#	get some reference geo info to put into the plot of 
#		tortoise coords.  this step is mostly just a sanity check
mojave.outlines <- map("county",c("california","nevada","arizona"),xlim=c(-120,-113),ylim=c(31,37))
mojave.outlines.sp <- map2SpatialLines(mojave.outlines,proj4string=CRS("+proj=longlat"))
mojave.outlines.sp2 <- spTransform(mojave.outlines.sp,cropped_lat@crs)

png(file="road_layer.png",res=150,width=7*150,height=7*150)
	plot(resampled_cropped_road,main="Distance from a Road")
	lines(mojave.outlines.sp2)
	points(proj.latlong.coords,pch=8)
dev.off()

################################
#	learning how to mask
################################

#	here, I'm figuring out how to mask a raster
#		by a given polygon.  that is, I want to take a raster
#		and a polygon that fits over some part of that raster
#		and I want to make all the frid cells in that raster that 
#		are outside of the polygon take the value 'NA'

#	first, I get the outline of the USA
#		and make it into a spatial polygon
#		in the same coordinate reference system as the rasters
usa.outlines <- map("usa")
IDs <- sapply(strsplit(usa.outlines$names, ":"), function(x) x[1])
usa.outlines.sp <- map2SpatialPolygons(usa.outlines,IDs=IDs,proj4string=CRS("+proj=longlat"))
usa.outlines.sp2 <- spTransform(usa.outlines.sp,cropped_lat@crs)

#	then, I mask the rasters by the polygon of the usa border, 
#		and check to see if it worked with some plotting
masked_resampled_cropped_road <- mask(resampled_cropped_road,usa.outlines.sp2,filename="masked_resampled_cropped_road")
masked_resampled_cropped_aspect <- mask(resampled_cropped_aspect,usa.outlines.sp2,filename="masked_resampled_cropped_aspect")
masked_cropped_precip <- mask(cropped_precip,usa.outlines.sp2,filename="masked_cropped_precip")
masked_cropped_lat <- mask(cropped_lat,usa.outlines.sp2,filename="masked_cropped_lat")
masked_cropped_long <- mask(cropped_long,usa.outlines.sp2,filename="masked_cropped_long")

png(file="masked_road_layer.png",res=150,width=7*150,height=7*150)
	plot(masked_resampled_cropped_road,main="Distance from a Road")
	lines(mojave.outlines.sp2)
	points(proj.latlong.coords,pch=8)
dev.off()

png(file="masked_aspect_layer.png",res=150,width=7*150,height=7*150)
	plot(masked_resampled_cropped_aspect,main="Distance from a Road")
	lines(mojave.outlines.sp2)
	points(proj.latlong.coords,pch=8)
dev.off()

################################
#	some earlier notes
################################

# this is a combination of stuff that Peter
#	and Gideon worked out about how rasters work
#	and how much memory they take

# check out what a raster is all about
str(road)

# the resolution of a raster determines the objects size, and
#	hence, the ability to read it into working memory and play around with it.
#	you can take a raster and decrease its resolution using the aggregate()
#	function, which takes as its arguments:
#			-the raster you want to zoom out of
#			-the factor by which you want to decrease resolution
#			-the function used to combine cell values of aggregated cells
#			-the file to output to
# caveat scriptor: this step can take a while and have high RAM requirements
system.time( aggregate(road,fact=10,fun=mean,na.rm=TRUE,filename="road_redux",overwrite=TRUE) )  # 45 secs, not much RAM usage?

# make a brick out of these?
system.time( rpa.stack <- stack( list(road,precip,aspect) ) )
# Oh, dear: Error in compareRaster(x) : different extent -- precip differs.
precip.cropped <- crop( precip, road, filename="precip_cropped", snap="in",overwrite=TRUE )
road.cropped <- crop( road, precip.cropped, filename="road_cropped", snap="in",overwrite=TRUE )
compareRaster(road.cropped,precip.cropped)  # still don't match??!?!?!!

precip.resample <- resample( precip.cropped, road.cropped, filename="precip_resample" )
compareRaster( precip.resample, road.cropped )

# ok, out of just two of them
system.time( rpa.stack <- stack( list(road,aspect) ) )
system.time( brick( rpa.stack, file="road_aspect_brick" ) )  # huh: two 2.2G files makes an 8G file.
rpa.brick <- brick("road_aspect_brick")

brick.cropped <- crop(rpa.brick, precip.cropped)

system.time( brick.redux <- aggregate(rpa.brick,fact=10,fun=mean,na.rm=TRUE,filename="brick_redux",overwrite=TRUE) )

precip.redux <- aggregate(precip.cropped,fact=10,fun=mean,na.rm=TRUE)
brick.redux.cropped <- crop( brick.redux, precip.cropped )
plot(brick.redux.cropped)
plot(precip.redux)
compareRaster(brick.redux.cropped,precip.redux)  # they DO match at the smaller resolution.

# NOT NECESSARY: remove other stuff from memory (but road isn't in memory, just a pointer to its file)
# rm(road) ; gc()

# now we can read in this aggregated road file, which, 
#	clocking in at only about 2300 x 2300, is much more manageable
road.redux <- raster("road_redux.grd")
plot(road.redux)


# get x,y coordinates for each cell
row.coords <- yFromRow(road.redux)
col.coords <- xFromCol(road.redux)

# now we read the raster into working memory as a matrix, 
#	dropping metadata along the way, and have a standard
#	R matrix object that we can do with as we wish
road.redux.matrix <- as.matrix(road.redux)

# what about the whole thing? 
if (FALSE) {
    # this uses 4.2G of RAM even though the files are only 2.2G. huh.
    #   Also note that it used significantly more ram while reading the thing in.
    system.time( road.matrix <- as.matrix(road) )  # took 36.760s
    sort( sapply(ls(),function(x){object.size(get(x))}))  # uses 4244121160
    rm(road.matrix); gc()
    # is that more than a just plain matrix?  nope.
    road.matrix <- matrix( 1.1, nrow=dim(road)[1], ncol=dim(road)[2] )  
    sort( sapply(ls(),function(x){object.size(get(x))}))  # uses 4244121160
    rm(temp.matrix); gc()
}

# then we can add x and y coordinates to the matrix
rownames(road.redux.matrix) <- row.coords
colnames(road.redux.matrix) <- col.coords

# and we can save that object to work with
save(file="road_redux_matrix.Robj",road.redux.matrix)

if(FALSE){
	# now, try to read in a geographic coordinate system
	lat <- raster("~/Downloads/evan_dt_gis/lat_30_utm_tif/lat_gcs_30.tif")
		aggregate(lat,fact=10,fun=mean,na.rm=TRUE,filename="lat_redux")
	long <- raster("~/Downloads/evan_dt_gis/lat_30_utm_tif/lon_gcs_30.tif")
		aggregate(long,fact=10,fun=mean,na.rm=TRUE,filename="long_redux")
}


# Memory profiling stuff:
sort( sapply(mget(ls()),object.size) )

# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}
# shorthand
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}



