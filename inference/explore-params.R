source("resistance-fns.R")
require(raster)
require(rgdal)
require(parallel)
numcores <- getcores()

config <- read.json.config("nussear-transforms/habitat-dem/config.json")
load("nussear-transforms/habitat-dem/setup.RData")

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
    for (k in 1:ncol(hts)) { ph(hts[,k],zlim=c(0,1e5)) }
    plot( p[1]+hts[locs,], pimat[,ref.inds], pch=20 ); abline(0,1)
    plot( p[1]+((hts[locs[ref.inds],]+t(hts[locs[ref.inds],]))/2), pimat[ref.inds,ref.inds], pch=20 ); abline(0,1)
    dev.off()
}

f(c(3e5,8,0,4,4,0,0,0))  # too long, no structure
f(c(3e5,8,2,4,4,0,4,0))  # a bit better
f(c(3e5,8,2,4,4,0,4,-5))  # hm, interesting. too long
f(c(3e5,8,2,4,4,0,8,-5))  # hm, interesting. too long
f(c(3e5,9,2,4,-4,-3,8,-5))  # hey, closer!
f(c(3e5,9,2,4,-4,-3,8,-7))  # hm, outliers
f(c(3e5,9,2,4,-4,-3,10,-7))  # got rid of most of outliers
f(c(3e5,8,3,4,-4,-3,10,-7))  # hm, ok-ish
f(c(3e5,8,3,4,-4,-3,10,-7))  # 
