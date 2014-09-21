################################################################
################################################################
#	go through all GIS layers to:
# -  make .pngs of them,
# -  count how many cells are NA in them
################################################################
################################################################

require(raster)
require(rgdal)
require(maps)
require(maptools)
require(sp)
rasterOptions(tmpdir=".")

if(!file.exists("masked_raster_pngs")){
	dir.create("masked_raster_pngs")
}
#	load the tortoise location data
#		made using "format_tortoise_location_data.R" in the git
load("tort.coords.rasterGCS.Robj")

county.outlines <- map("county",c("california","nevada","arizona"),xlim=c(-120,-113),ylim=c(31,37),plot=FALSE)
county.outlines <- map2SpatialLines(county.outlines,proj4string=CRS("+proj=longlat"))
county.outlines <- spTransform(county.outlines,tort.coords.rasterGCS@proj4string)

raster.tifs <- list.files()[grepl(".gri$",list.files())]

for(i in 1:length(raster.tifs)){
	current.raster <- raster(raster.tifs[i])
	raster.name <- strsplit(list.files()[grepl(".gri$",list.files())][i],split=".gri")[[1]]
	cat(raster.name,"\n")
	na.count <- freq(current.raster,value=NA)
	png(file=paste("masked_raster_pngs/",raster.name,".png",sep=""),res=200,width=10*200,height=10*200)
		plot(current.raster,main=raster.name)
		lines(county.outlines)
		points(tort.coords.rasterGCS@coords,pch=8,cex=0.7)
		legend(x="bottomleft",
				legend=c(paste("NA count: ",na.count,sep=""),
						 paste("NA prop: ", round(na.count/(dim(current.raster)[1] * dim(current.raster)[2]),3))))
	dev.off()
}

#	the code below was used to make the original raster images,
#		but now the document is being repurposed to make raster
#		images of the masked rasters, as well as to count the number 
#		of NA cells in each raster.  However, I hate to get rid 
#		of code, so it's in the graveyard now...

#	GRAVEYARD
if(FALSE){
load("county_tortoise_plotting_info.Robj")  # @gbradburd: how was this produced?

if (file.exists("geolayers/TIFF")) { setwd("geolayers/TIFF") }

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
	dev.off()},error= function(e){print(e);cat("problem with ",raster.name,"\n")})
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
}