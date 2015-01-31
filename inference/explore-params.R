source("../resistance-fns.R")
require(raster)
require(rgdal)
require(parallel)
numcores <- getcores()

config <- read.json.config("habitat-dem/config.json")
load("habitat-dem/setup.RData")

## look at where samples go
ph <- plot.ht.fn(nonmissing=nonmissing,layer=nalayer,sample.loc.file="../../tort_272_south_info/geog_coords.RData")

diag(pimat) <- NA
pimat[pimat<2e5] <- NA
ref.inds <- sample(which.nonoverlapping(neighborhoods),16)

f <- function (p) {
    r.inds <- sample(ref.inds,16)
    G@x <- update.G(p[-1])
    hts <- hitting.analytic(neighborhoods[r.inds],G,numcores=numcores)
    png(file="Rplots.png",width=24*144,height=12*144,pointsize=10,res=144)
    layout(matrix(1:18,nrow=3,byrow=TRUE))
    for (k in 1:ncol(hts)) { ph(hts[,k]) }
    plot( p[1]+hts[locs,], pimat[,r.inds], pch=20 ); abline(0,1)
    plot( p[1]+((hts[locs[r.inds],]+t(hts[locs[r.inds],]))/2), pimat[r.inds,r.inds], pch=20 ); abline(0,1)
    dev.off()
}

f(paramvec(config))

