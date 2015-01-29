#!/usr/bin/Rscript

require(raster)

nus <- raster("nussear.gri")

nus_gt_zero <- (nus>0)
writeRaster(nus_gt_zero,file="nussear_gt_zero.grd")

other.dir <- "../../expanded/64x"
for (x in c("agp_250","dem_30")) {
    r <- raster(file.path(other.dir,x))
    rp <- projectRaster(r,crs=CRS(proj4string(nus)))
    rs <- resample(rp,nus)
    stopifnot( compareRaster(nus,rs) )
    writeRaster(rs,file=paste(x,".grd",sep=''))
}
