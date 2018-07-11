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
wbd_adj <- lapply(seq_along(raw_wbd_adj), function (k) {
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

write.csv(tort_wbd, file="watershed_assignments.csv") # note this is modified again later on
for (k in seq_along(wbd_adj)) {
    write.csv(as.matrix(wbd_adj[[k]]), file=sprintf("WBD%d_adjacency.csv", 2*k))
}

if (FALSE) {
    # sanity checks
    wbd_num <- 4
    plot(wbd.list[[wbd_num]][nonempty_polys[[wbd_num]]])
    text(tort.coords, labels=tort_wbd[,wbd_num])
}

### Get a better WBD8-like map:

wbd_adj <- lapply(1:5, function (k) {
                read.csv(sprintf("WBD%d_adjacency.csv", 2*k), 
                         header=TRUE, check.names=FALSE, row.names=1) } )

wbd8 <- wbd.list[[4]][colnames(wbd_adj[[4]])]
wbd10 <- wbd.list[[5]]

merge_poly <- function (polys, merge_ids) {
    idarg <- names(polys)
    new_id <- paste(names(polys)[merge_ids], collapse=',')
    idarg[merge_ids] <- new_id
    # merged <- gUnaryUnion(polys[merge_ids], id=idarg)
    # SpatialPolygons(c(polys[-merge_ids]@polygons, merged@polygons))
    gUnaryUnion(polys, id=idarg)
}

mod_wbd8 <- wbd8
mod_wbd8 <- merge_poly(mod_wbd8, c(12, 15, 16, 19))
mod_wbd8 <- merge_poly(mod_wbd8, c(15, 16))
mod_wbd8 <- merge_poly(mod_wbd8, c(3,4))

plot(mod_wbd8, col=adjustcolor(rainbow(length(mod_wbd8)), 0.2))
points(tort.coords)
text(gCentroid(mod_wbd8, byid=TRUE), labels=seq_along(mod_wbd8))

split_poly <- function (polys, split) {
    split <- gUnaryUnion(split, id=rep(paste(names(split), collapse=","), length(split)))
    outside <- gDifference(polys, split, byid=c(TRUE,FALSE))
    SpatialPolygons(c(outside@polygons, split@polygons))
}

sub_wbd8 <- mod_wbd8
sub_wbd8 <- split_poly(sub_wbd8, wbd10[c(749, 750, 751)])
sub_wbd8 <- split_poly(sub_wbd8, wbd10[c(724, 677, 678)])
sub_wbd8 <- merge_poly(sub_wbd8, c(1,16))
sub_wbd8 <- merge_poly(sub_wbd8, c(1,8))
sub_wbd8 <- split_poly(sub_wbd8, wbd10[c(360)])
sub_wbd8 <- merge_poly(sub_wbd8, c(10,15))
sub_wbd8 <- split_poly(sub_wbd8, wbd10[c(330,334,566)])
sub_wbd8 <- merge_poly(sub_wbd8, c(9,15))
sub_wbd8 <- sub_wbd8[-8]
sub_wbd8 <- merge_poly(sub_wbd8, c(6,7))
sub_wbd8 <- split_poly(sub_wbd8, wbd10[c(424,409)])

tort_handpicked <- apply(gCovers(sub_wbd8, tort.coords, byid=TRUE), 1, which)
stopifnot(sort(unique(tort_handpicked)) == 1:length(sub_wbd8))
adj_handpicked <- adjacency(sub_wbd8)
nonempty_hp <- 1:length(sub_wbd8)
hp_adj <- sparseMatrix(i=rep(1:length(sub_wbd8), sapply(adj_handpicked, length)),
                        j=unlist(adj_handpicked),
                        x=1L, dims=c(length(sub_wbd8),length(sub_wbd8)))
rownames(hp_adj) <- colnames(hp_adj) <- names(sub_wbd8)

save(sub_wbd8, file="handpicked_WBD8.RData")
write.csv(as.matrix(hp_adj), file="handpicked_WBD8_adjacency.csv")
tort_wbd <- read.csv("watershed_assignments.csv", header=TRUE, row.names=1)
tort_wbd[,"handpicked_WBD8"] <- tort_handpicked
write.csv(tort_wbd, file="watershed_assignments.csv")

if (FALSE) {
    # sanity checks
    plot(sub_wbd8, col=adjustcolor(rainbow(length(sub_wbd8)), 0.2))
    points(tort.coords)
    text(gCentroid(sub_wbd8, byid=TRUE), labels=seq_along(sub_wbd8))

    plot(wbd10, col=adjustcolor(rainbow(length(wbd10)), 0.2))
    text(gCentroid(wbd10, byid=TRUE), labels=seq_along(wbd10))
}
