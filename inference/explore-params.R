
layer.prefix <- "../geolayers/multigrid/256x/crm_"
subdir <- "256x"
layer.file <- "dem-layer-list"
param.file <- "params-dem-layer-list.tsv"


source("resistance-fns.R")
require(raster)

require(parallel)
numcores <- getcores()

layer.names <- scan(layer.file,what="char") 

load(paste(subdir,"/",basename(layer.prefix),basename(layer.file),"-","setup.RData",sep=''))  # provides below, also 'pimat'
# load( paste(subdir,"/",basename(layer.prefix),"_",basename(layer.file),"_","G.RData",sep='') ) # provides "G"    "Gjj"    "update.G" "ndelta"   "ngamma"   "transfn"  "valfn"    "layers"
# load( paste( subdir, "/", basename(layer.prefix), basename(layer.file), "_neighborhoods.RData", sep='' ) ) # provides 'neighborhoods'
# load(paste(subdir,"/",basename(layer.prefix),"tortlocs.RData",sep='')) # provides 'locs'
# load( paste(subdir, "/", basename(layer.prefix),"_", basename(layer.file),"_nonmissing.RData",sep='') ) # provides nonmissing

# REMOVE MISSING INDIV
na.indiv <- which( is.na( locs ) )
locs <- locs[-na.indiv]
neighborhoods <- lapply(neighborhoods[-na.indiv],function (x) { x[!is.na(x)] })


##
# initial parameters?
#   in time t, 1D RW with rate r does t*r jumps, displacement has variance t*r
#   so time to move N grid sites away is sqrt(N)/r
#   so if hitting times are of order T, want r of order sqrt(N)/T

init.param.table <- read.table( param.file, header=TRUE )
init.params <- unlist( init.param.table[ match( subdir, init.param.table[,1] ), -1 ] )


ph <- plot.ht.fn(layer.prefix,"dem_30_m800_sq",nonmissing)
dem <- raster(paste(layer.prefix,"dem_30_m800_sq",sep=''))

dothese <- c(1,10,19,83)

newparams <- function (dothese) {
    G@x <- update.G(params)
    hts <- hitting.analytic( neighborhoods[dothese], G-diag(rowSums(G)), numcores=numcores )
    hts[hts<0] <- NA
    ymax <- 1.5*max(hts[locs,])
    for (k in seq_along(dothese)) { ph(pmin(ymax,hts[,k])); } # ph(pmin(6,log10(hts[,k])))  }
    invisible(hts)
}

layout(matrix(1:4,nrow=2))
params <- c(1e-4,-2,-9); hts <- newparams(c(10,83))
params <- c(.01,0,-1); hts <- newparams(c(10,83))
params <- c(.01,0,-3); hts <- newparams(c(10,83))
params <- c(.01,-.1,-3); hts <- newparams(c(10,83))

plot( exp((-2)*dem) )

# for talk
params <- c(.01,0,-1); hts1<- newparams(c(10,83))
params <- c(.01,0,-3); hts2 <- newparams(c(10,83))

ph <- plot.ht.fn(layer.prefix,"dem_30_m800_sq",nonmissing,par.args=list(mai=c(.7,.7,.1,.15),mar=c(3,3,0,0)+.1))

# see raster:::.rasterImagePlot

pdf(file="../../talks/images/example-hts.pdf", width=5, height=3.5, pointsize=10)
layout(matrix(1:4,nrow=2,byrow=TRUE))
for (k in 1:2) {
    ph(pmin(ymax,hts[,k]),legend.mar=0,legend=FALSE)
    ph(pmin(ymax,hts[,k]),legend.mar=0,legend=FALSE)
}
dev.off()


###
# find parameters that best fit the data?



