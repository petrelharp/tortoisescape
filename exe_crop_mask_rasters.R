################################################################
################################################################
#	Run crop_mask_raster on all specified rasters in directory
################################################################
################################################################

#	load the required packages
#	source the function code
#	and create the desired directories 
#		into which the output files will be saved
require(raster)
rasterOptions(tmpdir=".")
options(error=recover)
source("crop_mask_rasters.R")
dir.create("masked")
dir.create("10x")
dir.create("100x")

#	make a list of all the raster .tif files
#		that we want to read in, crop, resample,
#		mask, and aggregate by 10x and 100x

raster.tifs <- list.files()[grepl(".tif$",list.files())]

#	remove specific .tif files from the list

raster.tifs <- raster.tifs[-c(
								grep("m2_01tmax",raster.tifs),		# in GCS rather than UTM like the others
								grep("m2_01tmean",raster.tifs),		# in GCS rather than UTM like the others
								grep("dem_10",raster.tifs),			# didn't make a .png w/ raster_viz
								grep("surfarea_10",raster.tifs),	# didn't make a .png w/ raster_viz
								grep("surfratio_10",raster.tifs),	# didn't make a .png w/ raster_viz
								grep("avg_rough_30",raster.tifs)	# the smallest raster, which everything else is being cropped to
							)]

#	go through and crop, resample, mask, and aggregate the .tifs
#		in the specified list

for(i in 1:length(raster.tifs)){
	cat("starting on",raster.tifs[i],"\n")
	tryCatch({
		crop_mask_aggregate(raster2crop.file = raster.tifs[i],
							raster2crop.outputfilename = strsplit(raster.tifs[i],split=".tif")[[1]],
							reference_raster.file = "avg_rough_30.tif",
							masking_object.file = "usa_border_maskingObject.Robj")
			 },error=function(e){cat("problem with",raster.tifs[i],"\n")})
	removeTmpFiles()
}

