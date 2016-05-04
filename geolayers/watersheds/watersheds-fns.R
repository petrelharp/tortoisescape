
# from http://robinlovelace.net/r/2014/07/29/clipping-with-r.html
gClip <- function(shp, bb){
  if(class(bb) == "matrix") { 
      b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  } else {
      b_poly <- as(extent(bb), "SpatialPolygons")
  }
  gIntersection(shp, b_poly, byid = TRUE)
}

#' For each geometry in spgeom1, return the index of the geometry in spgeom2 that it overlaps with most.
hierarchy <- function ( spgeom1, spgeom2 ) {
    # this gives a list of length equal to that of spgeom2
    # of indices of spgeom1 that might overlap
    tree <- rgeos::gBinarySTRtreeQuery( spgeom2, spgeom1 )
    max.inds <- sapply( seq_along(spgeom1), function (k) {
            areas <- gArea( gIntersection( spgeom1[k], spgeom2[tree[[k]]], byid=c(FALSE,TRUE), drop_lower_td=TRUE ), byid=TRUE )
            match( names(areas)[which.max(areas)], names(spgeom2) )
        } )
    return(max.inds)
}

#' Return a list of color schemes that reflects a hierarchy
#' @param hier A list of numeric vectors, each giving the categories that the previous one are part of, so that every element of hier[[j]] is in 1...length(hier[[j-1]])
color_hierarchy <- function (hier, chroma=function(h){50}, luminance=function(h){50+50*runif(length(h))}, ...) {
    max.branches <- sapply( lapply( hier, table ), max ) + c(0,rep(1,length(hier)-1))
    # hues[[k]] gives a vector of locations that once modded to [0,1] correspond to the locations on the color circle of the corresponding nodes
    hues <- lums <- vector("list",length(hier))
    .lproj <- function (tt) { 50 + 50*abs( 2*((tt-50)/100 - floor((tt-50)/100+1/2))) }  # reflect into (50,100) at both endpoints
    for (j in seq_along(hier)) {
        hues[[j]] <- if (j>1) { hues[[j-1]][ hier[[j]] ] } else { rep(0,length(hier[[j]])) }
        for (k in unique(hier[[j]])) {
            dothese <- ( hier[[j]] == k )
            hues[[j]][dothese] <- hues[[j]][dothese] + ( seq_len(sum(dothese)) - sum(dothese)/2 )/prod(max.branches[1:j])
        }
    }
    for (j in seq_along(hier)) {
        lums[[j]] <- if (j>1) { lums[[j-1]][ hier[[j]] ] } else { rep(0,length(hier[[j]])) }
        for (k in unique(hier[[j]])) {
            dothese <- ( hier[[j]] == k )
            dlums <- floor( 15*sample(seq(-1,1,length.out=sum(dothese))) )
            lums[[j]][dothese] <- .lproj( lums[[j]][dothese] + dlums )
        }
    }
    # return( list( hues, lums ) )
    return( lapply( seq_along(hues), 
           function (j) { 
               hcl( 
                   h=360*hues[[j]], 
                   c=chroma(hues[[j]]), 
                   l=luminance(hues[[j]]), 
                   ... ) 
           } ) )
}

.expand <- function (x,fac) {
    mx <- median(x,na.rm=TRUE)
    mx + (1+fac)*(x-mx)
}
