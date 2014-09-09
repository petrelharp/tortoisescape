################################################################
################################################################
#	Run crop_mask_raster on all specified rasters in directory
################################################################
################################################################

require(raster)
rasterOptions(tmpdir=".")
options(error=recover)
source("crop_mask_rasters.R")
dir.create("masked")
dir.create("10x")
dir.create("100x")

raster.tifs <- list.files()[grepl(".tif$",list.files())]

raster.tifs <- raster.tifs[-c(
								grep("m2_01tmax",raster.tifs),
								grep("m2_01tmean",raster.tifs),
								grep("dem_10",raster.tifs),
								grep("surfarea_10",raster.tifs),
								grep("surfratio_10",raster.tifs),
								grep("avg_rough_30",raster.tifs)
							)]
							
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

