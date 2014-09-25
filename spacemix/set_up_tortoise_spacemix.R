#set up a spacemix run for the tortoises


require(raster)
require(rgdal)

#first, get the lat/long data down
load("tort.coords.rasterGCS.Robj")

tort.coords.GCS.spobj <- spTransform(tort.coords.rasterGCS,CRS("+proj=longlat"))
tort.coords <- tort.coords.GCS.spobj@coords

#second, get the sample covariance matrix

sample.covariance <- as.matrix(read.table("~/Desktop/Dropbox/tortoisescape/covmat/alleleCounts500kLoci-covmat.txt"))

#third, specify the mean sample sizes
mean.sample.sizes <- rep(1,nrow(sample.covariance))

save(tort.coords,sample.covariance,mean.sample.sizes,file="tort180_spacemix_dataset.Robj")