#!/usr/bin/Rscript

###
# Get hitting times with a landscape layer
#   e.g.
#     Rscript initial-hitting-times.R ../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_ annual_precip
# (note the space!)                                                                           ---->^

require(parallel)
numcores<-as.numeric(scan(pipe("cat /proc/cpuinfo | grep processor | tail -n 1 | awk '{print $3}'")))+1

source("resistance-fns.R")
require(raster)

if (!interactive()) {
    layer.prefix <- commandArgs(TRUE)[1]
    subdir <- commandArgs(TRUE)[2]
    layer.file <- commandArgs(TRUE)[3]
    param.file <- if (length(commandArgs(TRUE))>3) { commandArgs(TRUE)[4] } else { NULL }
} else {
    layer.prefix <- c("../geolayers/TIFF/500x/500x_")
    subdir <- "500x"
    layer.file <- "six-raster-list"
    param.file <- NULL
    # layer.names <- c("imperv_30", "agp_250", "m2_ann_precip", "avg_rough_30", "dem_30", "bdrock_ss2_st")
}
layer.names <- scan(layer.file,what="char") 

load(paste(subdir,"/",basename(layer.file),"-",basename(layer.prefix),"setup.RData",sep=''))

##
# initial parameters?
#   in time t, 1D RW with rate r does t*r jumps, displacement has variance t*r
#   so time to move N grid sites away is sqrt(N)/r
#   so if hitting times are of order T, want r of order sqrt(N)/T

gridwidth <- sqrt(dim(G)[1])  # roughly, N
ratescale <- sqrt(gridwidth)/mean(pimat)

if (is.null(param.file)) {
    init.params <- c( beta=ratescale, gamma=rep(1,length(layer.names)), delta=rep(1,length(layer.names) ) )
} else {
    init.params <- scan( param.file )
} 

G@x <- update.G(init.params)

###
# Conjugate gradient

dG <- rowSums(G)
cG <- colSums(G)
# objective function
H <- function (par,obs.ht,loc,locs,g.match=1) {
    a <- par[1]
    ht <- par[-1]
    ht[loc] <- 0
    z <- G%*%ht - dG*ht + 1
    z[loc] <- 0
    return( ( sum( z^2 ) + g.match * sum( (ht[locs] - (obs.ht-a) )^2 ) )/length(z) )
}
dH <- function (par,obs.ht,loc,locs,g.match=1) {
    # cG - G[loc,] is, except at [loc], 1^T ((G-diag(dG))[-loc,])
    a <- par[1]
    ht <- par[-1]
    z <- G%*%ht - dG*ht
    z[loc] <- 0
    z <- (G%*%z - dG*z) + (cG-G[loc,]) 
    z[locs] <- z[locs] + g.match*(ht[locs]-(obs.ht-a))
    z[loc] <- 0
    return( c( 2 * g.match * sum( ht[locs] - (obs.ht-a) ) / length(z) , 2 * as.vector(z) / length(z) ) )
}

# parscale <- rep( nrow(G) / exp( mean( log(dG), trim=.1, na.rm=TRUE ) ), nrow(G) )
parscale <- c( min(pimat), rep( mean(pimat), nrow(G) ) )
init.hts <- matrix(parscale,nrow=nrow(G)+1,ncol=length(locs))
init.hts[cbind(locs,seq_along(locs))] <- 0

# H( init.hts[,1], obs.ht=pimat[,1], loc=locs[1], locs=locs )
# dH( init.hts[,1], obs.ht=pimat[,1], loc=locs[1], locs=locs )

# k <- 28 
# test.ht <- optim( par=init.hts[,k], fn=H, gr=dH, obs.ht=pimat[,k], loc=locs[k], locs=locs, g.match=1/200,
#         method="L-BFGS-B", control=list( parscale=parscale, maxit=1000 ), lower=0, upper=Inf ) 
# ph(test.ht$par[-1]); with( environment(ph), points(tort.coords.rasterGCS[k]) )

optim.ht.list <- mclapply( seq_along(locs), function (k) {
            optim( par=init.hts[,k], fn=H, gr=dH, obs.ht=pimat[,k], loc=locs[k], locs=locs, 
                method="L-BFGS-B", control=list( parscale=parscale, maxit=1000 ), lower=0, upper=Inf ) 
        }, mc.cores=numcores )

convergences <- sapply(optim.ht.list,"[[","convergence")
unconverged <- which(convergences != 0)

for (k in 1:3) {
    if (length(unconverged)<length(locs)) {
        optim.ht.list[unconverged] <- mclapply( unconverged, function (k) {
                    newstart <- sample( setdiff(seq_along(locs),unconverged), 1 )
                    optim( par=optim.ht.list[[newstart]]$par, fn=H, gr=dH, obs.ht=pimat[,k], loc=locs[k], locs=locs, 
                        method="L-BFGS-B", control=list( parscale=parscale, maxit=1000 ), lower=0, upper=Inf ) 
                }, mc.cores=numcores )
        convergences <- sapply(optim.ht.list,"[[","convergence")
        unconverged <- which(convergences != 0)
    }
}

optim.hts <- sapply(optim.ht.list,"[[","par")

if (any(convergences!=0)) { warning("Some did not converge") }

save( layer.prefix, layer.names, subdir, optim.hts, convergences, init.params, ratescale, gridwidth, file=paste(subdir,"/",basename(layer.prefix),basename(layer.file),"-init-hts.RData",sep='') )


###
# testing
if (FALSE) {


    load(paste(subdir,"/",basename(layer.file),"-",basename(layer.prefix),"setup.RData",sep=''))
    load("../tort.coords.rasterGCS.Robj")
    load(paste(subdir,"/",basename(layer.prefix),"alllocs.RData",sep='')) # all.locs.dists
    load("torts-info.RData")
    ph <- plot.ht.fn(layer.prefix,"annual_precip",nonmissing)

    optim.hts[cbind(1+locs,seq_along(locs))] <- NA

    for (k in 1:ncol(optim.hts)) {
        layout(matrix(1:4,nrow=2))
        ph( optim.hts[-1,k], main=optim.hts[1,k] )
        points(tort.coords.rasterGCS[k])
        plot( all.locs.dists[,k], optim.hts[-1,k], pch=20, cex=.5 )
        points( all.locs.dists[locs,k], optim.hts[1+locs,k], col='red' )
        plot( tort.dists[k,], optim.hts[(1+locs),k], pch=20, cex=.5 )
        plot( pimat[,k]-optim.hts[1,k], optim.hts[(1+locs),k], pch=20, cex=.5 ); abline(0,1)
        if (is.null(locator(1))) { break }
    }


    ###
    # analytic

    tmp.pimat <- pimat-mean(diag(pimat))

    solve.hts <- interp.hitting( fullG, locs, tmp.pimat )
    solve.hts[cbind(locs,seq_along(locs))] <- 0

    plot( as.vector(tmp.pimat), as.vector(solve.hts[locs,]), col=1+(row(pimat)==col(pimat)) ); abline(0,1)

    for (k in seq_along(locs)[1:length(locs)]) {
        plot.ht( pmax(solve.hts[,k],0), hitting.layer, nonmissing )
        text( tort.coords.rasterGCS, labels=1:180 )
        points( tort.coords.rasterGCS[k+if(k>56){1}else{0}], pch="*", cex=4, col='red' )
        if (is.null(locator(1))) { break }
    }
}

