################################
#	getting raster data from tortoises locations
################################


#	this function extracts the raster value
#		from the matrix of locations provided
require(raster)

get.raster.data.from.locations <- function(raster,locations){
	values <- extract(raster, locations)
	return(values)
}


#	so, now we go through all the masked rasters
#		and apply this function to them.

load("tort.coords.rasterGCS.Robj")
raster.files <- list.files(pattern=".gri")
raster.names <- unlist(strsplit(raster.files,"crop_resampled_masked_"))[seq(2,2*length(raster.files),2)]
raster.names <- gsub(".gri","",raster.names)
raster.values.tort.locations.matrix <- matrix(0,
										nrow=nrow(tort.coords.rasterGCS@coords),
										ncol=length(raster.files))
row.names(raster.values.tort.locations.matrix) <- row.names(tort.coords.rasterGCS@coords)
colnames(raster.values.tort.locations.matrix) <- raster.names

for(i in 1:length(raster.files)){
	cat(i,"\t")
	tryCatch({
	tmp <- raster(raster.files[i])
		},error=function(e){cat("problem with",raster.files[i],"\n")})
		raster.values.tort.locations.matrix[,i] <- get.raster.data.from.locations(tmp,tort.coords.rasterGCS@coords)
}

save(raster.values.tort.locations.matrix,file="raster.values.tort.locations.matrix.Robj")
write.table(raster.values.tort.locations.matrix,file="raster.values.tort.locations.matrix.csv")