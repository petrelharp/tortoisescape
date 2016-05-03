################################################################
################################################################
#	Function to Crop, Mask, and Aggregate Rasters
################################################################
################################################################


# raster2crop.file is the filename of the raster you want to crop
# reference_raster.file is the filename of the raster to which you want to crop the raster2crop
# masking_object.file is the .Robj containing the polygon used to mask the raster
# raster2crop.outputfilename is the base name of all output files
#		-e.g., if raster2crop is in "road.tif", the raster2crop.outputfilename should be "road"
if(FALSE){
	#SAMPLE USAGE
	crop_mask_aggregate(raster2crop.file = "~/Downloads/tortTestGIS/road_30.tif",
						raster2crop.outputfilename = "road",
						reference_raster.file = "~/Downloads/tortTestGIS/lat_gcs_30.tif",
						masking_object.file = "usa_border_maskingObject.Robj",
						aggregate = 10)
}

crop_mask_aggregate <- function(raster2crop.file, raster2crop.outputfilename, reference_raster.file, masking_object.file){
	require(raster)
#		recover()
		# read in files as rasters
		raster2crop <- raster(raster2crop.file)
		reference_raster <- raster(reference_raster.file)

		# crop the focal raster to the reference raster,
		# and save a new, cropped raster file
		prefix <- "crop_"
		cat("cropping",raster2crop.file,"...\n")
		raster2crop <- crop(raster2crop,
							reference_raster,snap="in",
							filename=paste(prefix,raster2crop.outputfilename,sep=""),
							overwrite=TRUE)
		
		# if necessary, resample the focal raster to the 
		# reference raster
		if( xmax(raster2crop) != xmax(reference_raster) || 
			xmin(raster2crop) != xmin(reference_raster) ||
			ymax(raster2crop) != ymax(reference_raster) ||
			ymin(raster2crop) != ymin(reference_raster)){
			cat("resampling",raster2crop.file,"...\n")
			prefix <- paste(prefix,"resampled_",sep="")
			raster2crop  <- resample(raster2crop,
										reference_raster,
										filename=paste(prefix,raster2crop.outputfilename,sep=""),
										overwrite=TRUE)
		}
		
		# mask the cropped and/or resampled focal raster to 
		# the masking spatial polygon object
		cat("reading in masking object","\n")
		masking.object <- load(masking_object.file)
			load(masking_object.file)
		masking.object <- get(masking.object)
		prefix <- paste(prefix,"masked_",sep="")
		cat("masking",raster2crop.file,"...\n")
		raster2crop <- 	mask(raster2crop,
								masking.object,
								filename=paste("masked/",prefix,raster2crop.outputfilename,sep=""),
								overwrite=TRUE)

		# aggregate the focal raster 10X to
		# make a smaller raster object
		prefix10 <- paste(prefix,"aggregated_10x_",sep="")
		cat("aggregating 10x",raster2crop.file,"...\n")
			aggregate(raster2crop,
						fact=10,
						fun=mean,
						na.rm=TRUE,
						filename=paste("10x/",prefix10,raster2crop.outputfilename,sep=""),
						overwrite=TRUE)

		prefix100 <- paste(prefix,"aggregated_100x_",sep="")
		cat("aggregating 100x",raster2crop.file,"...\n")
			aggregate(raster2crop,
						fact=100,
						fun=mean,
						na.rm=TRUE,
						filename=paste("100x/",prefix100,raster2crop.outputfilename,sep=""),
						overwrite=TRUE)
		focal.raster.files <- list.files()[grep(raster2crop.outputfilename,list.files())]
		focal.raster.files <- focal.raster.files[grep("crop",focal.raster.files)]
		for(i in 1:length(focal.raster.files)){
			if(file.exists(focal.raster.files[i])){
				file.remove(focal.raster.files[i])
			}
		}
	return(0)
}

# GRAVEYARD
if(FALSE){
		# focal.raster.files <- focal.raster.files[-grep("masked",focal.raster.files)]
		#
		# if specified, the reference raster will be cropped 
		# back to the newly cropped focal raster.  this should
		# only be necessary if the reference raster is not
		# already cropped.
		if(!is.null(reciprocal_crop)){
			writeRaster(
				crop(raster2,
					extent(
						max(raster2@extent@xmin,
								cropped_raster1@extent@xmin),
						min(raster2@extent@xmax,
								cropped_raster1@extent@xmax),
						max(raster2@extent@ymin,
								cropped_raster1@extent@ymin),
						min(raster2@extent@ymax,
								cropped_raster1@extent@ymax))),
					filename=paste("crop",reference_raster.file,sep=""),overwrite=TRUE)
		reference_raster <- raster(paste("crop",reference_raster.file,sep=""))
		}
}