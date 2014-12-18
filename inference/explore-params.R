
layer.prefix <- "../geolayers/multigrid/256x/crm_"
subdir <- "256x"
layer.file <- "dem-layer-list"


source("resistance-fns.R")
require(raster)

require(parallel)
numcores <- getcores()

layer.names <- scan(layer.file,what="char") 

load(paste(subdir,"/",basename(layer.prefix),basename(layer.file),"-","setup.RData",sep=''))  # provides below, also 'pimat'
# load( paste(subdir,"/",basename(layer.prefix),"_",basename(layer.file),"_","G.RData",sep='') ) # provides "G"    "Gjj"    "update.G" "ndelta"   "ngamma"   "transfn"  "valfn"    "layers"
# load( paste( subdir, "/", basename(layer.prefix), basename(layer.file), "_neighborhoods.RData", sep='' ) ) # provides 'neighborhoods'
# load(paste(subdir,"/",basename(layer.prefix),basename(layer.file),"_tortlocs.RData",sep='')) # provides 'locs'
# load( paste(subdir, "/", basename(layer.prefix),"_", basename(layer.file),"_nonmissing.RData",sep='') ) # provides nonmissing

## LOGISTIC FUNCTION
transfn <- function (x) { 1/(1+exp(-x)) }

# REMOVE MISSING INDIV
na.indiv <- which( is.na( locs ) )
locs <- locs[setdiff(seq_along(locs),na.indiv)]
neighborhoods <- lapply(neighborhoods[setdiff(seq_along(locs),na.indiv)],function (x) { x[!is.na(x)] })


## END SETUP

ph <- plot.ht.fn(layer.prefix,"dem_30_m800_sq",nonmissing)
# ph( layers[,1] )

dothese <- c(1,10,19,83)

newparams <- function (params,dothese,do.layout=TRUE) {
    # params are parameters: for n layers,
    #  params[1] is beta, overall multiplicative constant
    #  params[2:(n+1)] is gamma, weights on the layers that give the stationary distribution
    #  params[(n+2):(2*n+1)] is delta, the weights on the layers that give the jump rates
    #
    # dothese is a vector of indices of tortoise locations to compute hitting times of
    #
    # Produces 2*n+4 plots.
    if (do.layout) {
        nplots <- 2*length(dothese)+4
        layout(matrix(1:(floor(sqrt(nplots))*ceiling(sqrt(nplots))),nrow=floor(sqrt(nplots)),byrow=TRUE))
    }
    G@x <- update.G(params)
    hts <- hitting.analytic( neighborhoods[dothese], G-diag(rowSums(G)), numcores=numcores )
    hts[hts<0] <- NA
    # 'locs' are indices of tortoise locations
    ymax <- 1.5*max(hts[locs,])  # 1.5 times maximum hitting time to another tortoise location
    # plot with maximum value at ymax
    for (k in seq_along(dothese)) { ph(pmin(ymax,hts[,k]), main=paste("hitting time to ", dothese[k]) ) }
    plot( hts[locs,], pimat[,dothese], col=col(pimat[,dothese]), xlab='hitting time', ylab='divergence' )
    abline(0,1,untf=TRUE)
    plot( hts[locs,], pimat[,dothese], col=col(pimat[,dothese]), xlab='hitting time', ylab='divergence', log='xy' )
    ph( valfn( params[1 + (1:ngamma)] ), main="stationary distribution" )
    ph( valfn( params[1 + ngamma + (1:ndelta)] ), main="relative jump rate" )
    invisible(hts)
}

hts <- newparams(c(1e-4,-2,-9),c(10,83))


hts <- newparams(c(.01,0,-1),c(10,83))

hts <- newparams(c(.01,0,-3),c(10,83))

hts <- newparams(c(.01,-.1,-3),c(10,83))


# what does a transformed layer look like?
dem <- raster(paste(layer.prefix,"dem_30_m800_sq",sep=''))
plot( exp((-2)*dem) )

