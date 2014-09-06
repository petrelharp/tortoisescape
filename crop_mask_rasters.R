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
	crop_mask_aggregate(raster2crop.file = "road.tif",
						raster2crop.outputfilename = "road",
						reference_raster.file = "crop_lat.gri",
						masking_object.file = "usa_border_maskingObject.Robj",
						aggregate = 10)
}

crop_mask_aggregate <- function(raster2crop.file, raster2crop.outputfilename, reference_raster.file, masking_object.file, aggregate=NULL){
	require(raster)
		# read in files as rasters
		raster2crop <- raster(raster2crop.file)
		reference_raster <- raster(reference_raster.file)
		
		# crop the focal raster to the reference raster,
		# and save a new, cropped raster file
		prefix <- "crop_"
		writeRaster(
			crop(raster2crop,
				extent(
					max(reference_raster@extent@xmin,
							raster2crop@extent@xmin),
					min(reference_raster@extent@xmax,
							raster2crop@extent@xmax),
					max(reference_raster@extent@ymin,
							raster2crop@extent@ymin),
					min(reference_raster@extent@ymax,
							raster2crop@extent@ymax))),
				filename=paste(prefix,raster2crop.outputfilename,sep=""),overwrite=TRUE)
		raster2crop <- raster(paste(prefix,raster2crop.outputfilename,sep=""))
		
		# if necessary, resample the focal raster to the 
		# reference raster
		if( xmax(raster2crop) != xmax(reference_raster) || 
			xmin(raster2crop) != xmin(reference_raster) ||
			ymax(raster2crop) != ymax(reference_raster) ||
			ymin(raster2crop) != ymin(reference_raster)){
			cat("resampling...","\n")
			prefix <- paste(prefix,"resampled_",sep="")
			resample(raster2crop,
						reference_raster,
						filename=paste(prefix,raster2crop.outputfilename,sep=""),
						overwrite=TRUE)
			raster2crop <- raster(paste(prefix,raster2crop.outputfilename,sep=""))
		}
		
		# mask the cropped and/or resampled focal raster to 
		# the masking spatial polygon object
		masking.object <- load(masking_object.file)
			load(masking_object.file)
		masking.object <- get(masking.object)
		prefix <- paste(prefix,"masked_",sep="")
		writeRaster(
			mask(raster2crop,
					masking.object,
					filename=paste(prefix,raster2crop.outputfilename,sep="")),
			overwrite=TRUE)
		raster2crop <- raster(paste(prefix,raster2crop.outputfilename,sep=""))
		
		# if specified, aggregate the focal raster to
		# make a smaller raster object
		if(!is.null(aggregate)){
			prefix <- paste(prefix,"aggregated_",sep="")
			writeRaster(
				aggregate(raster2crop,
							fact=aggregate,
							fun=mean,
							na.rm=TRUE,
							filename=paste(prefix,raster2crop.outputfilename,sep=""),
							overwrite=TRUE))
		}
	return(0)
}

# GRAVEYARD
if(FALSE){
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