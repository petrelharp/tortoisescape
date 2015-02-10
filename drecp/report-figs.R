#!/usr/bin/Rscript

require(raster)
require(rgdal)

dem <- raster("../geolayers/nussear/habitat-model/dem_30.grd")
nus <- raster("../geolayers/nussear/habitat-model/nussear.grd")
crew <- raster("../geolayers/nussear/habitat-model/crew_NA_1km.gri")
bb <- extent( 3.4e5, 8.7e5, 3.62e6, 4.17e6 )

tort.loc.obj <- load("../tort_272_info/geog_coords.RData")
assign( "sample.locs", spTransform( get(tort.loc.obj), CRSobj=CRS(proj4string(nus)) ) )
ref.pt.obj <- load("../geolayers/reference_points/all_ref_points.RData")
assign( "ref.points", spTransform( get(ref.pt.obj), CRSobj=CRS(proj4string(nus)) ) )
isolated.refs <- c( 2, 16, 44, 50, 136, 142, 153, 174 )

png(width=3*288, height=3*288, res=288, file="reference-points.png", pointsize=10)
#pdf(width=2.5, height=2.5, file="reference-points.pdf", pointsize=10)
par(mar=c(1,1,3,1)+.1)
plot( crop(mask(nus,crew),bb), legend=FALSE,xaxt='n',yaxt='n',main='reference points')
points(ref.points[-isolated.refs], pch=20, cex=0.25)
dev.off()
