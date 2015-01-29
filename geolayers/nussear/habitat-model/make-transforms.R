#!/usr/bin/Rscript

require(raster)

nus <- raster("nussear.gri")

nus_gt_zero <- (nus>0)
writeRaster(nus_gt_zero,file="nussear_gt_zero.grd")
nus_gt_one <- (nus>0.1)
writeRaster(nus_gt_one,file="nussear_gt_one.grd")
nus_gt_two <- (nus>0.2)
writeRaster(nus_gt_two,file="nussear_gt_two.grd")
nus_gt_three <- (nus>0.3)
writeRaster(nus_gt_three,file="nussear_gt_three.grd")

other.dir <- "../../expanded/64x"
for (x in c("agp_250","dem_30","lat_utm_30","lon_utm_30")) {
    r <- raster(file.path(other.dir,x))
    rp <- projectRaster(r,crs=CRS(proj4string(nus)))
    rs <- resample(rp,nus)
    stopifnot( compareRaster(nus,rs) )
    writeRaster(rs,file=paste(x,".grd",sep=''),overwrite=TRUE)
}
