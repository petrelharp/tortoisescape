source("resistance-fns.R")
require(raster)
require(rgdal)
require(parallel)
numcores <- getcores()

config <- read.json.config("habitat/config.json")
load("habitat/setup.RData")

## look at where samples go
ph <- plot.ht.fn(nonmissing=nonmissing,layer=nalayer)

diag(pimat) <- NA
pimat[pimat<2e5] <- NA
ref.inds <- c(10,23,42,179)
ref.inds <- c(21,24,3,105,213,97,114,235,35,222,227,15,12,89,155,179)

f <- function (p) {
    G@x <- update.G(p[-1])
    hts <- hitting.analytic(neighborhoods[ref.inds],G,numcores=numcores)
    png(file="Rplots.png",width=24*144,height=12*144,pointsize=10,res=144)
    layout(matrix(1:18,nrow=3,byrow=TRUE))
    for (k in 1:ncol(hts)) { ph(hts[,k]) }
    plot( p[1]+hts[locs,], pimat[,ref.inds], pch=20 ); abline(0,1)
    plot( p[1]+((hts[locs[ref.inds],]+t(hts[locs[ref.inds],]))/2), pimat[ref.inds,ref.inds], pch=20 ); abline(0,1)
    dev.off()
}

f(c(3e5,8,0,8,0,0))  # no structure
f(c(3e5,5,0,8,0,3))  # looks promising but too long
f(c(3e5,5,0,8,-1,2)) # looks promising but too long
f(c(3e5,5,6,8,-1,2)) # right scale! but too much structure
f(c(2.9e5,5,6,8,-1,1)) # right scale! needs more structure
f(c(2.9e5,5,6,8,-3,5)) # similar; has outliers
f(c(2.9e5,5,4,10,-2,5)) # too long, not enough structure
f(c(2.9e5,7,4,10,-2,5)) # too long, not enough structure
