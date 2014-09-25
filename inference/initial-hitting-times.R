#!/usr/bin/Rscript

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
# now follow interp.hitting in resistance-fns.R

Pmat <- sparseMatrix( i=seq_along(locs), j=locs, x=1, dims=c(length(locs),nrow(G)) )
PtP <- crossprod(Pmat)
dG <- rowSums(G)
interp.hts <- sapply( seq_along(locs), function (kk) {
            Gk <- G[-locs[kk],]
            diag(Gk) <- (-1)*dG[-locs[kk]]
            bvec <- crossprod(Pmat,pimat[,kk]) - crossprod( Gk, rep(1.0,nrow(G)-1) )
            as.numeric( solve( PtP+crossprod(Gk), bvec ) )
} )

tmp <- interp.hitting( G - diag(rowSums(G)), locs, pimat )

################
# previous attempt
if (FALSE) {

# get some initial values for the iterative solver
load(paste(basename(layer.prefix),"alllocs.RData",sep='')) # provides all.locs.dists
all.locs.dists <- all.locs.dists[,-na.indiv]

ht.lms <- lapply( seq_along(locs), function (kk) {
        lm( pimat[kk,-kk] ~ dz, data.frame(dz=all.locs.dists[locs[-kk],kk]) )
    } )
init.hts <- sapply( seq_along(ht.lms), function (kk) {
        predict( ht.lms[[kk]], newdata=data.frame(dz=all.locs.dists[,kk]), )
    } )

jacobi.hts.list <- vector(mode="list",length=nreps)
jacobi.hts.list[[1]] <- init.hts 

for (k in 2:nreps) {
    jacobi.hts.list[[k]] <- hitting.jacobi(locs,G,jacobi.hts.list[[k-1]],kmax=kmax)
}

save( jacobi.hts.list, init.params, ratescale, gridwidth, file=paste(basename(layer.prefix),layer.name,"-init-hts.RData",sep='') )

###
# look at results
if (FALSE) {

    kk <- 1
    x <- sapply( jacobi.hts.list, function (y) y[,kk] )
    x[locs[kk],] <- 0
    Gx <- G%*%x - rowSums(G)*x
    Gx[locs[kk],] <- 0

    matplot(Gx,pch=20,cex=0.25)

}
}
