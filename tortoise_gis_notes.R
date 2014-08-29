################################################################
################################################################
#	NOTES ON WORKING WITH GIS LAYERS IN R	
################################################################
################################################################
require(raster)

# this put the temp stuff in my temp directory, and filled up that partition:
rasterOptions(tmpdir=".")

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
road <- raster("tortTestGIS/road_30.tif")
precip <- raster("tortTestGIS/annual_precip.tif")
aspect <- raster("tortTestGIS/aspect_30.tif")

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
