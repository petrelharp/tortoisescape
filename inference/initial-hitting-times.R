#!/usr/bin/Rscript

require(parallel)

###
# Get hitting times with a landscape layer
#   e.g.
#     Rscript initial-hitting-times.R ../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_ annual_precip
# (note the space!)                                                                           ---->^

nreps <- 3
kmax <- 1e4

source("resistance-fns.R")
require(raster)

if (!interactive()) {
    layer.prefix <- commandArgs(TRUE)[1]
    layer.name <- commandArgs(TRUE)[2]
    subdir <- if (length(commandArgs(TRUE))>2) { commandArgs(TRUE)[3] } else { "." }
} else {
    layer.prefix <- c("../geolayers/TIFF/500x/500x_")
    layer.name <- "annual_precip"
    subdir <- "."
}

# get precomputed G
load(paste(subdir,"/",basename(layer.prefix),"G.RData",sep=''))
load(paste(subdir,"/",basename(layer.prefix),"nonmissing.RData",sep=''))
Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )

###
# layer whatnot

layers <- cbind( scale( values( raster(paste(layer.prefix,layer.name,sep='')) )[nonmissing] ) )
stopifnot(nrow(layers)==nrow(G))

# tortoise locations
load(paste(subdir,"/",basename(layer.prefix),"tortlocs.RData",sep=''))
nind <- length(locs)
na.indiv <- which( is.na( locs ) )
locs <- locs[-na.indiv]

# pairwise divergence values
pimat.vals <- scan("../pairwisePi/alleleCounts_1millionloci.pwp") # has UPPER with diagonal
pimat <- numeric(nind^2)
dim(pimat) <- c(nind,nind)
pimat[upper.tri(pimat,diag=TRUE)] <- pimat.vals
pimat[lower.tri(pimat,diag=FALSE)] <- t(pimat)[lower.tri(pimat,diag=FALSE)]
pimat <- pimat[-na.indiv,-na.indiv]

# scale to actual pairwise divergence, and then by 1/mutation rate
pimat <- pimat * .018 * 1e8

##
# initial parameters?
#   in time t, 1D RW with rate r does t*r jumps, displacement has variance t*r
#   so time to move N grid sites away is sqrt(N)/r
#   so if hitting times are of order T, want r of order sqrt(N)/T

gridwidth <- sqrt(dim(G)[1])  # roughly, N
ratescale <- sqrt(gridwidth)/mean(pimat)

init.params <- c( beta=ratescale, gamma=1, delta=1 )

G@x <- update.G(init.params)

###
# Conjugate gradient

dG <- rowSums(G)
cG <- colSums(G)
# objective function
H <- function (ht,obs.ht,loc,locs,g.match=1) {
    ht[loc] <- 0
    z <- G%*%ht - dG*ht + 1
    z[loc] <- 0
    return( ( sum( z^2 ) + g.match * sum( (ht[locs] - obs.ht)^2 ) )/length(z) )
}
dH <- function (ht,obs.ht,loc,locs,g.match=1) {
    # cG - G[loc,] is, except at [loc], 1^T ((G-diag(dG))[-loc,])
    z <- G%*%ht - dG*ht
    z[loc] <- 0
    z <- (G%*%z - dG*z) + (cG-G[loc,]) 
    z[locs] <- z[locs] + g.match*(ht[locs]-obs.ht)
    z[loc] <- 0
    return( 2 * as.vector(z) / length(z) )
}

# parscale <- rep( nrow(G) / exp( mean( log(dG), trim=.1, na.rm=TRUE ) ), nrow(G) )
parscale <- rep( mean(pimat), nrow(G) )
init.hts <- matrix(parscale,nrow=nrow(G),ncol=length(locs))
init.hts[cbind(locs,seq_along(locs))] <- 0

optim.ht.list <- mclapply( seq_along(locs), function (k) {
            optim( par=init.hts[,k], fn=H, gr=dH, obs.ht=pimat[,k], loc=locs[k], locs=locs, method="L-BFGS-B", control=list( parscale=parscale ), lower=0, upper=Inf ) 
        } )
optim.hts <- sapply(optim.ht.list,"[[","par")

if (any(sapply(optim.ht.list,"[[","convergence")!=0)) { warning("Some did not converge") }

save( optim.hts, init.params, ratescale, gridwidth, file=paste(subdir,"/",basename(layer.prefix),layer.name,"-init-hts.RData",sep='') )


###
# testing
if (FALSE) {

###
# analytic

tmp.pimat <- pimat-mean(diag(pimat))

solve.hts <- interp.hitting( fullG, locs, tmp.pimat )
solve.hts[cbind(locs,seq_along(locs))] <- 0

plot( as.vector(tmp.pimat), as.vector(solve.hts[locs,]), col=1+(row(pimat)==col(pimat)) ); abline(0,1)

for (k in seq_along(locs)[1:length(locs)]) {
    plot.ht( (solve.hts[,k]), hitting.layer, nonmissing )
    text( tort.coords.rasterGCS, labels=1:180 )
    points( tort.coords.rasterGCS[k+if(k>56){1}else{0}], pch="*", cex=4, col='red' )
    if (is.null(locator(1))) { break }
}

for (k in seq_along(locs)[1:length(locs)]) {
    plot.ht( pmax(solve.hts[,k],0), hitting.layer, nonmissing )
    text( tort.coords.rasterGCS, labels=1:180 )
    points( tort.coords.rasterGCS[k+if(k>56){1}else{0}], pch="*", cex=4, col='red' )
    if (is.null(locator(1))) { break }
}

    load("../tort.coords.rasterGCS.Robj")

    fullG <- G
    diag(fullG) <- (-1)*rowSums(G)
    true.hts <- hitting.analytic(locs,fullG)
    hlayer <- raster(paste(layer.prefix,layer.name,sep='')) 
    values(hlayer)[-nonmissing] <- NA # NOTE '-' NOT '!'
    k <- 1
    values(hlayer)[nonmissing] <- true.hts[,k]
    plot(hlayer)

    solve.hts <- interp.hitting( fullG, locs, true.hts[locs,] )
    solve.hts[cbind(locs,seq_along(locs))] <- 0

    range(solve.hts)

    k <- 1
    H(solve.hts[,k], obs.ht=solve.hts[locs,k], loc=k, locs=locs )
    dH(solve.hts[,k], obs.ht=solve.hts[locs,k], loc=k, locs=locs )

}

