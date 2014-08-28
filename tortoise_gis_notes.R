################################################################
################################################################
#	NOTES ON WORKING WITH GIS LAYERS IN R	
################################################################
################################################################
require(raster)
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
road <- raster("~/Downloads/tortTestGIS/road_30.tif")

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
aggregate(road,fact=10,fun=mean,na.rm=TRUE,filename="road_redux")

# OPTIONAL: remove other stuff from memory
rm(road) ; gc()

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