source("../resistance-fns.R") 
require(raster)
require(rgdal)
require(parallel)
numcores <- getcores()

config <- read.json.config("dem/config.json")
load("dem/setup.RData")

## look at where samples go
ph <- plot.ht.fn(nonmissing=nonmissing,layer=nalayer,sample.loc.file="../../tort_272_info/geog_coords.RData")

diag(pimat) <- NA
pimat[pimat<2e5] <- NA
ref.inds <- which.nonoverlapping(neighborhoods)

f <- function (p) {
    G@x <- update.G(p[-1])
    hts <- hitting.analytic(neighborhoods[ref.inds],G,numcores=numcores)
    png(file="Rplots.png",width=20*144,height=12*144,pointsize=10,res=144)
    layout(matrix(1:15,nrow=3,byrow=TRUE))
    for (k in 1:ncol(hts)) { ph(hts[,k],zlim=c(0,1e5)) }
    plot( p[1]+hts[locs,], pimat[,ref.inds], pch=20 ); abline(0,1)
    plot( p[1]+((hts[locs[ref.inds],]+t(hts[locs[ref.inds],]))/2), pimat[ref.inds,ref.inds], pch=20 ); abline(0,1)
    dev.off()
}

f(paramvec(config))
