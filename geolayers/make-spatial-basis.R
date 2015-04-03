require(raster) 
require(Matrix)

base.dir <- gsub("tortoisescape/.*","tortoisescape",getwd())

source(file.path(base.dir,"inference","resistance-fns.R"))
source(file.path(base.dir,"inference","input-output-fns.R"))

nalayer <- raster("nussear/habitat-model/mask_crew_dem_2K_sea_habitat")
nonmissing <- which(!is.na(values(nalayer)))
ph <- plot.ht.fn(layer.prefix="nussear/habitat-model/",nonmissing=nonmissing,

ref.pt.file <- file.path(base.dir,"geolayers/reference_points/all_ref_points.RData")
ref.pt.obj  <- load(ref.pt.file)
ref.points <- get(ref.pt.obj)

ref.dists <- do.call( cbind, lapply( seq_along(ref.points), function (k) {
            values( distanceFromPoints(nalayer,ref.points[k]) )[!is.na(values(nalayer))] 
        } ) )


ph( ref.dists[,1] )

