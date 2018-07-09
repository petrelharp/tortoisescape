# Find which tortoises are in which watershed.
# Run this in e.g. tort_272_info

library(rgdal)
library(raster)
library(sp)
library(rgeos)
library(spdep)
library(colorspace)
library(maps)
source("../visualization/map-utils.R",chdir=TRUE)
watershed_dir <- "../geolayers/watersheds"
source(file.path(watershed_dir, "watersheds-fns.R"), chdir=TRUE)

tort.coords.obj <- load("geog_coords.RData")
tort.coords <- get(tort.coords.obj)

dem <- get_elev()
conts <- get_contours(dem)
shade <- get_shading(dem)
counties <- get_counties(dem)

# note this is not doing WBD 12, which needs more RAM
wbd.list <- lapply( 2*(1:5), function (wbd_num) {
            gClip( 
                  spTransform( 
                               readOGR(file.path(watershed_dir, 
                                                 "WBD/WBD_National_931.gdb"),paste0("WBDHU",wbd_num)), 
                           CRS(proj4string(dem)) ), 
                     dem )
    } )

tort_wbd_mats <- lapply(wbd.list, function (wbd) {
                            gCovers(wbd, tort.coords, byid=TRUE) } )

stopifnot(all(unlist(lapply(tort_wbd_mat, rowSums)) == 1))

raw_tort_wbd <- do.call(cbind, lapply(tort_wbd_mats, function (x) apply(x, 1, which)))
colnames(raw_tort_wbd) <- paste0("WBD", 2*(1:5))


###
# get pairwise adjacency matrices

# Finds the self-adjacency matrix of the polygons
# defined to be any pair of polygons within distance eps of each other
adjacency <- function (polys) {
    spdep::poly2nb(polys, foundInBox=rgeos::gUnarySTRtreeQuery(polys))
}
# in the adjacency list, '0' means no neighbors (the channel islands)
raw_wbd_adj <- lapply(wbd.list, adjacency)

###
# Renumber polygons to refer to only ones that have tortoises

nonempty_polys <- apply(raw_tort_wbd, 2, unique)
wbd_adj <- lapply(seq_along(wbd_adj), function (k) {
                      adjlist <- lapply(raw_wbd_adj[[k]], setdiff, 0)
                      stopifnot(length(adjlist) == length(raw_wbd_adj[[k]]))
                      n <- max(unlist(adjlist))
                      adj <- sparseMatrix(i=rep(1:n, sapply(adjlist, length)),
                                          j=unlist(adjlist),
                                          x=1L, dims=c(n,n))
                      return(adj[nonempty_polys[[k]], nonempty_polys[[k]]])
                } )
for (k in seq_along(wbd_adj)) {
    rownames(wbd_adj[[k]]) <- colnames(wbd_adj[[k]]) <- names(wbd.list[[k]])[nonempty_polys[[k]]]
}

tort_wbd <- sapply(1:ncol(raw_tort_wbd), function (k) {
                       match(raw_tort_wbd[,k], nonempty_polys[[k]]) } )
dimnames(tort_wbd) <- dimnames(raw_tort_wbd)

write.csv(tort_wbd, file="watershed_assignments.csv")
for (k in seq_along(wbd_adj)) {
    write.csv(as.matrix(wbd_adj[[k]]), file=sprintf("WBD%d_adjacency.csv", 2*k))
}

if (FALSE) {
    # sanity checks
    wbd_num <- 4
    plot(wbd.list[[wbd_num]][nonempty_polys[[wbd_num]]])
    text(tort.coords, labels=tort_wbd[,wbd_num])
}
