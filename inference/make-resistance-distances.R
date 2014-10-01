#!/usr/bin/Rscript

###
# Get hitting times with a landscape layer
#   e.g.
#     Rscript initial-hitting-times.R ../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_ annual_precip
# (note the space!)                                                                           ---->^

source("resistance-fns.R")
require(raster)

require(parallel)
numcores<-as.numeric(scan(pipe("cat /proc/cpuinfo | grep processor | tail -n 1 | awk '{print $3}'")))+1

if (!interactive()) {
    layer.prefix <- commandArgs(TRUE)[1]
    subdir <- commandArgs(TRUE)[2]
    layer.file <- commandArgs(TRUE)[3]
    param.file <- if (length(commandArgs(TRUE))>3) { commandArgs(TRUE)[4] } else { NULL }
} else {
    # layer.prefix <- "../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_"
    # subdir <- "100x"
    # layer.file <- "../inference/six-raster-list"
    # param.file <- "500x/six-raster-list-hitting-times.tsv"

    layer.prefix <- "../geolayers/TIFF/500x/500x_"
    subdir <- "500x"
    layer.file <- "../inference/six-raster-list"
    param.file <- "simple-init-params-six-raster-list.tsv"
}
method <- "analytic"
layer.names <- scan(layer.file,what="char") 

load( paste(subdir,"/",basename(layer.prefix),"G.RData",sep='') ) # provides "G"        "update.G" "ndelta"   "ngamma"   "transfn"  "valfn"    "layers"
Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )

load(paste(subdir,"/",basename(layer.prefix),"tortlocs.RData",sep='')) # provides 'locs'
na.indiv <- which( is.na( locs ) )
locs <- locs[-na.indiv]

##
# initial parameters?
#   in time t, 1D RW with rate r does t*r jumps, displacement has variance t*r
#   so time to move N grid sites away is sqrt(N)/r
#   so if hitting times are of order T, want r of order sqrt(N)/T

init.param.table <- read.table( param.file, header=TRUE )
init.params <- as.vector( init.param.table[ match( subdir, init.param.table[,1] ), -1 ] )

G@x <- update.G(init.params)

if (method=="analytic") {

    hts <- hitting.analytic( locs, G-diag(rowSums(G)), numcores=numcores )

} else if (method=="CG") {

    dG <- rowSums(G)
    cG <- colSums(G)
    # objective function
    H <- function (ht,loc) {
        # ( (G-D)ht + 1 )^T S ( (G-D)ht + 1)
        # where S = I except S[loc,loc]=0
        ht[loc] <- 0
        z <- G%*%ht - dG*ht + 1
        z[loc] <- 0
        return( ( sum( z^2 ) )/length(z) )
    }
    dH <- function (ht,loc) {
        # 2 (G-D)^T S ( (G-D) ht + 1 )
        ht[loc] <- 0
        z <- G%*%ht - dG*ht + 1
        z[loc] <- 0
        z <- (crossprod(G,z) - dG*z)
        z[loc] <- 0
        return( 2 * as.vector(z) / length(z) )
    }

    # parscale <- rep( nrow(G) / exp( mean( log(dG), trim=.1, na.rm=TRUE ) ), nrow(G) )
    parscale <- rep( 1, nrow(G) )
    init.hts <- matrix(parscale,nrow=nrow(G),ncol=length(locs))
    init.hts[cbind(locs,seq_along(locs))] <- 0

    optim.ht.list <- mclapply( seq_along(locs), function (k) {
                optim( par=init.hts[,k], fn=H, gr=dH, loc=locs[k], 
                    method="L-BFGS-B", control=list( parscale=parscale, maxit=1000 ), lower=0, upper=Inf ) 
            }, mc.cores=numcores )

    convergences <- sapply(optim.ht.list,"[[","convergence")
    unconverged <- which(convergences != 0)

    if (FALSE) {
        # check gradient
        k <- 10
        loc <- locs[k]
        ht <- init.hts[,k]

        eps <- 1e-5 * runif(length(ht))
        c( H(ht,loc=loc), H(ht+eps,loc=loc)-H(ht,loc=loc), sum(eps*dH(ht,loc=loc)) )
    }
}

colnames( hts ) <- locs
write.table( hts, file=paste( subdir, "/", basename(layer.file), "-hitting-times.tsv", sep=''), row.names=FALSE )

if (FALSE) {

    load( paste(subdir, "/", basename(layer.prefix),"nonmissing.RData",sep='') ) # provides nonmissing
    ph <- plot.ht.fn(layer.prefix,"annual_precip",nonmissing)

    layout( matrix(1:6,nrow=2) )
    for (k in 1:ncol(hts)) { ph( hts[,k] ); if ((k%%6==0) && is.null(locator(1))) break }

}
