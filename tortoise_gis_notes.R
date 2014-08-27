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
road <- raster("road_30.tif")

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

# now we read the raster into working memory as a matrix, 
#	dropping metadata along the way, and have a standard
#	R matrix object that we can do with as we wish
road.redux.mat <- as.matrix(road.redux)

	