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
} else {
    layer.prefix <- c("../geolayers/TIFF/500x/500x_")
    layer.name <- "annual_precip"
}

# get precomputed G
load(paste(basename(layer.prefix),"G.RData",sep=''))
load(paste(basename(layer.prefix),"nonmissing.RData",sep=''))
Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )

###
# layer whatnot

layers <- cbind( scale( values( raster(paste(layer.prefix,layer.name,sep='')) )[nonmissing] ) )
stopifnot(nrow(layers)==nrow(G))

# tortoise locations
load(paste(basename(layer.prefix),"tortlocs.RData",sep=''))
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
H <- function (ht,obs.ht,loc,locs) {
    ht[loc] <- 0
    z <- G%*%ht - dG*ht
    z[loc] <- 0
    return( ( sum( z^2 ) + sum( (ht[locs] - obs.ht)^2 ) )/length(z) )
}
dH <- function (ht,obs.ht,loc,locs) {
    # cG - G[loc,] is, except at [loc], 1^T ((G-diag(dG))[-loc,])
    z <- G%*%ht - dG*ht
    z[loc] <- 0
    z <- (G%*%z - dG*z) + (cG-G[loc,]) 
    z[locs] <- z[locs] + (ht[locs]-obs.ht)
    z[loc] <- 0
    return( 2 * as.vector(z) / length(z) )
}

optim.hts <- optim( par=solve.hts[,k], fn=H, gr=dH, obs.ht=solve.hts[locs,k], loc=k, locs=locs, method="CG", control=list(parscale=rep(mean(solve.hts),nrow(G)),maxit=1000) )

optim.hts <- mclapply( seq_along(locs), function (k) {
            optim( par=init.hts[,k], fn=H, gr=dH, obs.ht=pimat[,k], locs=locs, method="CG", control=list( parscale=parscale ) ) 
        } )


###
# testing
if (FALSE) {

    fullG <- G
    diag(fullG) <- (-1)*rowSums(G)

    true.hts <- hitting.analytic(locs,fullG)
    hitting.layer <- raster(paste(layer.prefix,layer.name,sep='')) 
    values(hitting.layer)[-nonmissing] <- NA # NOTE '-' NOT '!'
    k <- 1
    values(hitting.layer)[nonmissing] <- true.hts[,k]
    plot(hitting.layer)

    solve.hts <- interp.hitting( fullG, locs, true.hts[locs,] )
    solve.hts[cbind(locs,seq_along(locs))] <- 0

    range(solve.hts)

    k <- 1
    H(solve.hts[,k], obs.ht=solve.hts[locs,k], loc=k, locs=locs )
    dH(solve.hts[,k], obs.ht=solve.hts[locs,k], loc=k, locs=locs )

}

save( jacobi.hts.list, init.params, ratescale, gridwidth, file=paste(basename(layer.prefix),layer.name,"-init-hts.RData",sep='') )
