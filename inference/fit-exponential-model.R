#!/usr/bin/Rscript

source("resistance-fns.R")
require(raster)

layer.prefix <- c("../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_")

# get precomputed G
load(paste(basename(layer.prefix),"G.RData",sep=''))
load(paste(basename(layer.prefix),"nonmissing.RData",sep=''))
Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )

###
# layer whatnot

layer.names <- c("annual_precip","barren_30","eastness_30","lat_gcs_30","lon_gcs_30")
layers <- sapply(layer.names, function (ll) {
            rast <- raster(paste(layer.prefix,ll,sep=''))
            # note this is ROW-ORDERED
            # so to plot do:  dim(x) <- dim(rast)[2:1]; image(x)
            vrast <- scale( values(rast)[nonmissing] )
            return(vrast)
        } )
stopifnot(nrow(layers)==nrow(G))

rm(nonmissing)

# tortoise locations
load(paste(basename(layer.prefix),"tortlocs.RData",sep=''))
na.indiv <- which( is.na( locs ) )
locs <- locs[-na.indiv]
nind <- length(locs)


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

init.params <- c( beta=ratescale, gamma=rep(.01,length(layer.names)), delta=rep(.01,length(layer.names)) )

G@x <- update.G(init.params)


# get some initial values for the iterative solver
load(paste(basename(layer.prefix),"alllocs.RData",sep='')) # provides all.locs.dists
all.locs.dists <- all.locs.dists[,-na.indiv]

ht.lms <- lapply( seq_along(locs), function (kk) {
        lm( pimat[kk,-kk] ~ dz, data.frame(dz=all.locs.dists[locs[-kk],kk]) )
    } )
init.hts <- sapply( seq_along(ht.lms), function (kk) {
        predict( ht.lms[[kk]], newdata=data.frame(dz=all.locs.dists[,kk]), )
    } )


system.time( init.hts <- hitting.jacobi(locs,G,init.hts,tol=.01,kmax=10) )

jacobi.hts <- hitting.jacobi(locs,G,init.hts,tol=.01,kmax=10)

jacobi.hts <- hitting.jacobi(locs,G,100*jacobi.hts,tol=.01,kmax=10)

jacobi.hts <- hitting.jacobi(locs,G,10*jacobi.hts,tol=.01,kmax=10)

jacobi.hts <- hitting.jacobi(locs,G,2*jacobi.hts,tol=.01,kmax=10)

if (FALSE) {
    load(paste(basename(layer.prefix),"nonmissing.RData",sep=''))
    rast <- raster(paste(layer.prefix,layer.names[1],sep=''))
    tmp <- matrix(NA,nrow=dim(rast)[2],ncol=dim(rast)[1])
    rm(rast)

    kk <- 1
    
    jhts <- jacobi.hts[,kk]

    tmp[nonmissing] <- jhts
    image(tmp)

    jhts <- hitting.jacobi( locs[kk], G, cbind(jhts) )

    plot(jhts,jacobi.hts[,kk]); abline(0,1)

    plot( G[-locs[kk],-locs[kk]] %*% jhts[-locs[kk]] - rowSums(G[-locs[kk],-locs[kk]])*jhts[-locs[kk]] )

}


