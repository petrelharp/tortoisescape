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


nus <- raster("../nussear/habitat-model/nussear.grd")
tort.loc.obj <- load("../../tort_272_info/geog_coords.RData")
assign( "sample.locs", spTransform( get(tort.loc.obj), CRSobj=CRS(proj4string(nus)) ) )

thresh <- 0.5
nsamps <- 1e3
samps <- xyFromCell( nus, which(values(nus)>thresh)[sample.int( sum(values(nus)>thresh,na.rm=TRUE), nsamps )] )
plot(nus>thresh)
points(samps,pch=20)
nonmissing <- which(!is.na(values(nus)))
neighborhoods <- get.neighborhoods( 1000, samps, nonmissing, nus )
nonoverlapping <- which.nonoverlapping(neighborhoods)

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
