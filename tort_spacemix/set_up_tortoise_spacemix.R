#set up a spacemix run for the tortoises


require(raster)
require(rgdal)

#	For the 270 tortoise dataset
#	first, get the lat/long data down
load("tort_272_info/geog_coords.RData")

tort.coords.GCS.spobj <- spTransform(tort.coords.272.rasterGCS,CRS("+proj=longlat"))
tort.coords <- tort.coords.GCS.spobj@coords

#	and drop the last 2 torts, as they are dupes
tort.coords <- tort.coords[-c(271:272),]

#	second, get the sample covariance matrix
#	drop the last 6 rows/columns, as they are from dupe torts
#		(the last 4 prepped w/ a different PCR protocol)

sample.covariance <- as.matrix(read.table("276torts.covar"))[-c(271:276),-c(271:276)]

#	name the rows and columns
row.names(sample.covariance) <- row.names(tort.coords)
colnames(sample.covariance) <- row.names(tort.coords)


#third, specify the mean sample sizes
mean.sample.sizes <- rep(100,nrow(sample.covariance))

tort270_spacemix_dataset <- list("tort.coords" = tort.coords,
								"sample.covariance" = sample.covariance,
								"mean.sample.sizes" = mean.sample.sizes)
save(tort270_spacemix_dataset,file="tort_spacemix/tort270_spacemix_dataset.Robj")





#For the 180 tortoise dataset
#first, get the lat/long data down
load("tort_180_info/tort.coords.rasterGCS.Robj")

tort.coords.GCS.spobj <- spTransform(tort.coords.rasterGCS,CRS("+proj=longlat"))
tort.coords <- tort.coords.GCS.spobj@coords

#second, get the sample covariance matrix

sample.covariance <- as.matrix(read.table("~/Desktop/Dropbox/tortoisescape/covmat/alleleCounts500kLoci-covmat.txt"))

#third, specify the mean sample sizes
mean.sample.sizes <- rep(100,nrow(sample.covariance))

save(tort.coords,sample.covariance,mean.sample.sizes,file="tort180_spacemix_dataset.Robj")
