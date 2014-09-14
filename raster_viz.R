################################################################
################################################################
#	go through all GIS layers to:
# -  make .pngs of them,
# - 
################################################################
################################################################

require(raster)
require(rgdal)
require(maps)
require(maptools)
require(sp)
rasterOptions(tmpdir=".")
load("county_tortoise_plotting_info.Robj")

raster.tifs <- list.files()[grepl(".tif$",list.files())]
raster.extents <- vector(mode="list",length=length(raster.tifs))

for(i in 1:length(raster.tifs)){
	current.raster <- raster(raster.tifs[i])
	raster.name <- strsplit(list.files()[grepl(".tif$",list.files())][i],split=".tif")[[1]]
	cat(raster.name,"\n")
	raster.extents[[i]] <- extent(current.raster)
		names(raster.extents)[i] <- raster.name
	tryCatch({
	png(file=paste(raster.name,".png",sep=""),res=120,width=10*120,height=10*120)
		par('plt' = c(0.1171429, 0.9400000, 0.1457143, 0.8828571))
		plot(current.raster,main=raster.name)
		lines(county_lines)
		points(tortoise_locations,pch=8)
	dev.off()},error= function(e){cat("problem with ",raster.name,"\n")})
	gc()
}

save(raster.extents,file="raster.extents.Robj")

raster.extents <- raster.extents[-c(
									which(
										names(raster.extents) == "m2_01tmax"  |
										names(raster.extents) =="m2_01tmean"  |
										names(raster.extents) =="dem_10" 	  |
										names(raster.extents) =="surfarea_10" |
										names(raster.extents) =="surfratio_10")
									)]


east.minz <- numeric(length(raster.extents))
south.minz <- numeric(length(raster.extents))
west.minz <- numeric(length(raster.extents))
north.minz <- numeric(length(raster.extents))
for(i in 1:length(raster.extents)){
cat(i,"\t")
	east.minz[i] <- raster.extents[[i]]@xmax
	south.minz[i] <- raster.extents[[i]]@ymin
	west.minz[i] <- raster.extents[[i]]@xmin
	north.minz[i]  <- raster.extents[[i]]@ymax
}

names(raster.extents)[which.min(east.minz)]
names(raster.extents)[which.max(west.minz)]
names(raster.extents)[which.min(north.minz)]
names(raster.extents)[which.max(south.minz)]
