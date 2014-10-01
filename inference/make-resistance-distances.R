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
    param.file <- commandArgs(TRUE)[4] 
    method <- commandArgs(TRUE)[5] 
    prev.ht <- if (length(commandArgs(TRUE))>5) { commandArgs(TRUE)[6] } else { NULL } 
} else {
    # layer.prefix <- "../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_"
    # subdir <- "100x"
    # layer.file <- "../inference/six-raster-list"
    # param.file <- "simple-init-params-six-raster-list.tsv"
    # method <- "CG"
    # prev.ht <- "100x/500x-aggregated-hitting-times.tsv"

    layer.prefix <- "../geolayers/TIFF/500x/500x_"
    subdir <- "500x"
    layer.file <- "../inference/six-raster-list"
    param.file <- "simple-init-params-six-raster-list.tsv"
    method <- "analytic"
    prev.ht <- NULL
}

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
init.params <- unlist( init.param.table[ match( subdir, init.param.table[,1] ), -1 ] )

G@x <- update.G(init.params)

if (method=="analytic") {

    hts <- hitting.analytic( locs, G-diag(rowSums(G)), numcores=numcores )

} else if (method=="CG") {

    init.hts <- as.matrix( read.table(prev.ht,header=TRUE) )

    stopifnot( nrow(init.hts) == nrow(G) )

    # infer better parameters

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
    parscale <- rep( mean(init.hts), nrow(init.hts) )

    optim.ht.list <- mclapply( seq_along(locs), function (loc.ind) {
                optim( par=init.hts[,loc.ind], fn=H, gr=dH, loc=locs[loc.ind], 
                    method="L-BFGS-B", control=list( parscale=parscale, maxit=1000 ), lower=0, upper=Inf ) 
            }, mc.cores=numcores )

    convergences <- sapply(optim.ht.list,"[[","convergence")
    unconverged <- which(convergences != 0)

    if (FALSE) {
        # check gradient
        loc.ind <- 10
        loc <- locs[loc.ind]
        ht <- init.hts[,loc.ind]

        eps <- 1e-5 * runif(length(ht))
        c( H(ht,loc=loc), H(ht+eps,loc=loc)-H(ht,loc=loc), sum(eps*dH(ht,loc=loc)) )

        # look at convergence
        load( paste(subdir, "/", basename(layer.prefix),"nonmissing.RData",sep='') ) # provides nonmissing
        ph <- plot.ht.fn(layer.prefix,"annual_precip",nonmissing)

        loc.ind <- 10
        oht.list <- vector( mode='list', length=6 )
        oht.list[[1]] <- list( par=init.hts[,loc.ind] )
        for (k in 2:10) { oht.list[[k]] <- optim( par=oht.list[[k-1]]$par, fn=H, gr=dH, loc=locs[loc.ind], method="L-BFGS-B", control=list( parscale=parscale, maxit=10000 ), lower=0, upper=Inf )  }

        layout( matrix(1:6,nrow=2) )
        for (k in 2:length(oht.list)) { 
            ph( oht.list[[k]]$par )
            ph( oht.list[[k]]$par - oht.list[[k-1]]$par )
            ph( dH( oht.list[[k-1]]$par, loc=locs[loc.ind] ) )
            if (is.null(locator(1))) break
        }

        layout( matrix(1:6,nrow=2) )
        ph( dH( oht.list[[length(oht.list)]]$par, loc=locs[loc.ind] ) )
        for (k in floor(seq(1,length(oht.list)-1,length.out=5))) {
            ph( oht.list[[length(oht.list)]]$par - oht.list[[k]]$par, main=k )
        }


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
