require(raster) 
require(Matrix)
require(sp)

base.dir <- gsub("tortoisescape/.*","tortoisescape",getwd())

source(file.path(base.dir,"inference","resistance-fns.R"))
source(file.path(base.dir,"inference","input-output-fns.R"))

nalayer <- raster("nussear/habitat-model/mask_crew_dem_2K_sea_habitat")
nonmissing <- which(!is.na(values(nalayer)))
ph <- plot.ht.fn(layer.prefix="nussear/habitat-model/",nonmissing=nonmissing)

ref.pt.file <- file.path(base.dir,"geolayers/reference_points/all_ref_points.RData")
ref.pt.obj  <- load(ref.pt.file)
ref.points <- get(ref.pt.obj)

ref.dists <- do.call( cbind, lapply( seq_along(ref.points), function (k) {
            values( distanceFromPoints(nalayer,ref.points[k]) )[!is.na(values(nalayer))] 
        } ) )


ph( ref.dists[,1] )


#######
# basis for local population size?

# choose something for the truth to check we can get close
popsize <- 1e4 * (0.05 + 2/(1+exp(-ref.dists[,1]*10/max(ref.dists))))
ph(popsize)

# subset
rpts <- ref.points[sample.int(length(ref.points),20)]
# pairwise distances
rpt.dists <- spDists(rpts)

midpoints <- vector(choose(length(rpts),2),mode="list")
n <- 1
for (k in seq_along(rpts)[-1]) {
    for (j in 1:(k-1)) {
        midpoints[[n]] <- (1/2)*( coordinates(rpts[k]) + coordinates(rpts[j]) )
        n <- n+1
    }
}

