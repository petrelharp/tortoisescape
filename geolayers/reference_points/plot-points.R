require(raster)
require(maptools)


nus <- raster("../nussear/habitat-model/nussear.grd")
plot(nus)

n.pts <- readShapePoints("north_regular.shp")
s.pts <- readShapePoints("south_regular.shp")
a.pts <- readShapePoints("all_regular.shp")
