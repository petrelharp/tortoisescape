require(raster)
require(rgdal)
require(maptools)
require(parallel)

get.neighborhoods <- function ( ndist, locations, nonmissing, layer, numcores=detectCores(), na.rm=TRUE ) {
    # locations should be either a SpatialPoints object or a 2-column matrix of coordinates
    this.lapply <- if ( numcores>1 && "parallel" %in% .packages()) { function (...) { mclapply( ..., mc.cores=numcores ) } } else { lapply }
    if ( class(locations)=="SpatialPoints" ) { locations <- coordinates(locations) }  # this is the first thing distanceFromPoints does anyhow
    if (is.null(dim(locations))) { locations <- matrix(locations,ncol=2) }
    neighborhoods <- this.lapply( 1:NROW(locations) , function (k) {
            d_tort <- distanceFromPoints( layer, locations[k,] ) 
            match( Which( d_tort <= max(ndist,minValue(d_tort)), cells=TRUE, na.rm=TRUE ), nonmissing )
        } )
    if (na.rm) { neighborhoods <- lapply(neighborhoods,function (x) { x[!is.na(x)] }) }
    return(neighborhoods)
}

which.nonoverlapping <- function (neighborhoods) {
    # find a set of neighborhoods that are mutually nonoverlapping
    perm <- sample(length(neighborhoods))
    goodones <- rep(FALSE,length(perm))
    goodones[1] <- TRUE
    for (k in seq_along(perm)[-1]) {
        goodones[k] <- ( 0 == length( intersect( neighborhoods[[k]], unlist(neighborhoods[goodones]) ) ) )
    }
    return( which(goodones) )
}

remove.clumps <- function (layer) {
    cl <- clump(layer,directions=4)
    big.clump <- which.max( table( values(cl) ) )
    layer[cl!=big.clump] <- NA
    return(layer)
}


nus <- raster("../nussear/habitat-model/nussear.grd")
tort.loc.obj <- load("../../tort_272_info/geog_coords.RData")
assign( "sample.locs", spTransform( get(tort.loc.obj), CRSobj=CRS(proj4string(nus)) ) )

masked <- raster("../nussear/habitat-model/mask_crew_dem_2K_sea.gri")
north <- raster("../nussear/habitat-model/mask_crew_dem_2K_sea_north.gri")
south <- raster("../nussear/habitat-model/mask_crew_dem_2K_sea_south.gri")

# the whole range
thresh <- 0.2
nus.thresh <- remove.clumps(mask(masked,nus>thresh,maskvalue=FALSE))
nsamps <- 5e3
neighbor.dist <- 10000
samps <- xyFromCell( nus.thresh, which(!is.na(values(nus.thresh)))[sample.int( sum(!is.na(values(nus))), nsamps )] )

nonmissing <- which(!is.na(values(nus)))
neighborhoods <- get.neighborhoods( neighbor.dist, samps, nonmissing, nus )
nonoverlapping <- which.nonoverlapping(neighborhoods)
ref.locs <- na.omit( samps[nonoverlapping,] )

ref.points <- SpatialPoints(coords=ref.locs,proj4string=CRS(proj4string(nus.thresh)),bbox=bbox(nus.thresh))

save(ref.points,file="all_ref_points.RData")

pdf(file="all-ref-points.pdf",width=5,height=5,pointsize=10)
plot(nus.thresh)
points(samps,pch=20,cex=0.5)
points(samps[nonoverlapping,],col='red')
dev.off()


# Jannet's:

n.pts <- readShapePoints("north_regular.shp")
s.pts <- readShapePoints("south_regular.shp")
a.pts <- readShapePoints("all_regular.shp")
a.jt.pts <- readShapePoints("all_int_jt.shp")  
n.mj.pts <- readShapePoints("north_int_mj.shp")  
s.jt.mj.pts <- readShapePoints("south_int_jt_mj.shp")

n.vals <- nus[cellFromXY( nus, n.pts )]
s.vals <- nus[cellFromXY( nus, s.pts )]
a.vals <- nus[cellFromXY( nus, a.pts )]
a.jt.vals <- nus[cellFromXY( nus, a.jt.pts )]
n.mj.vals <- nus[cellFromXY( nus, n.mj.pts )]
s.jt.mj.vals <- nus[cellFromXY( nus, s.jt.mj.pts )]
range(n.vals,s.vals,a.vals,a.jt.vals,n.mj.vals,s.jt.mj.vals)

pdf(file="point-locs.pdf",width=10,height=5,pointsize=10)
layout(t(1:2))
plot(nus)
points(sample.locs,pch=".")
points(n.pts,col='blue',pch=20)
points(s.pts,col='red',pch=20)
points(a.pts,pch=20)
plot(nus)
points(sample.locs,pch=".")
points(a.jt.pts,pch=20)
points(n.mj.pts,col='blue',pch=20)
points(s.jt.mj.pts,col='red',pch=20)
dev.off()
