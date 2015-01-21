#!/usr/bin/Rscript

# Make some simple layers on a small grid for testing purposes
set.seed(2357)

require(raster)
require(Matrix)
source("../resistance-fns.R")
require(parallel); numcores <- getcores()

# biglayer <- raster(nrows=2^10,ncols=2^10,xmn=0,ymn=0,xmx=100,ymx=100)
# values(biglayer) <- rnorm(length(biglayer))
# nlayers <- 2
# fact <- 2^5
# layer.list <- lapply( 1:nlayers, function (k) {
#         values(biglayer) <- rnorm(length(biglayer))
#         aggregate( biglayer, fact=2^5, fun=mean, na.rm=FALSE )
#     } )

nlayers <- 2
data(volcano)
v <- raster( (volcano-mean(volcano))/sd(volcano), xmn=2667400, xmx=2668010, ymn=6478700, ymx=6479570, crs="+init=epsg:27200")
v2 <- raster( (volcano^2-mean(volcano^2))/sd(volcano^2), xmn=2667400, xmx=2668010, ymn=6478700, ymx=6479570, crs="+init=epsg:27200")
layer.list <- list( v, v2 )

# include some missing values
navals <- sample(length(values(layer.list[[1]])),10)
for (k in 1:nlayers) { values(layer.list[[k]])[navals] <- NA }
layer <- layer.list[[1]]


# sampling locations
nlocs <- 40
SP.locs <- SpatialPoints( cbind( 
        sample(seq(bbox(layer)[1,1],bbox(layer)[1,2],length.out=100),nlocs),
        sample(seq(bbox(layer)[2,1],bbox(layer)[2,2],length.out=100),nlocs)
        ) )

nonmissing <- which(!is.na(values(layer)))

orig.locs <- cellFromXY( layer, SP.locs )
locs <- match(orig.locs,nonmissing)

ndist <- 10

neighborhoods <- get.neighborhoods( ndist, SP.locs, nonmissing, layer, numcores )
boundaries <- get.boundaries( neighborhoods, nonmissing, layer, numcores )
nonoverlapping <- which.nonoverlapping( neighborhoods )

layers <- sapply( layer.list, function (ll) { values(ll)[nonmissing] } )

G <- make.G( layer=layer, nonmissing=nonmissing )
Gjj <- p.to.j(G@p)

transfn <- function (x) { 1/(1+exp(-x)) }
valfn <- function (gamma) { ( rowSums( layers * gamma[col(layers)], na.rm=TRUE ) ) }
ndelta <- ngamma <- nlayers
update.G <- function(params) {
    beta <- params[1]
    gamma <- params[1+(1:ngamma)]
    delta <- params[1+ngamma+(1:ndelta)]
    return( exp(beta) * transfn(valfn(gamma))[G@i+1L] * transfn( valfn(delta)[G@i+1L] + valfn(delta)[Gjj] ) )
    # return( exp(beta) * transfn(valfn(gamma))[G@i+1L] * transfn( abs(valfn(delta)[G@i+1L] - valfn(delta)[Gjj]) ) )
}


true.params <- c(3*runif(1),8*(runif(2*nlayers)-0.5))
names(true.params) <- c("beta",paste("gamma",1:nlayers,sep=''),paste("delta",1:nlayers,sep=''))

G@x <- update.G(true.params)

tort.coords.rasterGCS <- SP.locs
save(tort.coords.rasterGCS, file="rasters/tort.coords.rasterGCS.Robj")

dir.create("rasters", showWarnings=FALSE)
for (k in seq_along(layer.list)) {
    writeRaster(layer.list[[k]], file=paste("rasters/layer_", letters[k], sep=''), overwrite=TRUE)
}

setup.list <- c( "G", "Gjj", "layers", "locs", "neighborhoods", "ndelta", "ngamma", "nonmissing", "transfn", "update.G", "valfn" )
save( list=setup.list, file="simple-setup.RData" )
