#!/usr/bin/Rscript

source("resistance-fns.R")

# width of square grid
n <- 100
nmap <- matrix(1:n^2,nrow=n)
all.locs <- data.frame( i=as.vector(row(nmap)), j=as.vector(col(nmap)) )

# random landscape
n.layers <- 2
AA <- replicate( n.layers, grid.generator(n), simplify=FALSE )
aa <- rep(1/n.layers,length=n.layers)
G <- aa[1] * AA[[1]]
for (k in 2:n.layers) {
    G <- G + aa[k] * AA[[k]]
}

# sampling locations
nsamps <- 100
locs.ij <- cbind( i=sample.int(n,nsamps)-1L, j=sample.int(n,nsamps)-1L )
locs <- ij.to.k(locs.ij,n,n)

# analytical mean hitting times
true.hts <- hitting.analytic(locs,G)  # warning, this could take a while (10s for n=100 and nsamps=20)

kk <- sample.int(nsamps,1)
plot( row(nmap), col(nmap), pch=20, col=colorize(true.hts[,kk]) )
points( locs.ij[kk,1], locs.ij[kk,2], pch="*", cex=2 )
points( locs.ij[,1], locs.ij[,2], pch="+", cex=1 )


# Observed hitting times
pairwise.hts <- true.hts[locs,]

# look at IBD:
pairwise.dists <- grid.dist( locs.ij )
plot( as.vector( pairwise.dists ), as.vector( pairwise.hts ), col=rainbow(n*1.2)[row(pairwise.dists)] )

# interpolate these. (note: gets NAs, which is fine... see 'all.zeros' below.)
interp.loess <- lapply( seq_along(locs), function (kk) {
            z <- log1p( pairwise.hts[,kk] )
            loess( z ~ i * j, data=locs.ij )
        } )
interp.hts <- sapply( interp.loess, function (interp) {
            expm1( predict( interp, newdata=all.locs ) )
        } )

# only a few percent error:
mean( (interp.hts-true.hts)/(true.hts+5), na.rm=TRUE )

# how's this do?
layout(matrix(1:40,nrow=8))
par(mar=c(0,0,0,0))
colorbreaks <- colorize( c(unlist(true.hts),unlist(interp.hts)), return.breaks=TRUE )
invisible( lapply(1:20, function (kk) { 
                    plot( j~i, col=colorize(true.hts[,kk],breaks=colorbreaks), data=all.locs, xlab='', ylab='', xaxt='n', yaxt='n' )
                    plot( j~i, col=colorize(interp.surf.hts[,kk],breaks=colorbreaks), data=all.locs, xlab='', ylab='', xaxt='n', yaxt='n' )
                } ) )


# now try to infer the coefficients:
#  want a,b to minimize 
#     | (a*A1+b*A2) %*% interp.hts + ones |^2

# these ones are the 'self hitting times' that we don't care about
zeros <- locs + (0:(length(locs)-1))*nrow(interp.hts)

# check this is right
ones <- (-1) * G %*% true.hts
range(ones[-zeros]) # yup
ones[zeros]  # huh, these are big though

# and how is interp doing
i.ones <- (-1) * G %*% interp.hts
hist(i.ones[-zeros],breaks=100) # whoa, not so great
i.ones[zeros]  # huh, these are big though

##
# ok, estimation
# think about scaling these things down

all.zeros <- unique( sort( c( zeros, which( is.na(i.ones) ) ) ) )

estimate.aa <- function (hts) {
    # ok, math now
    # here are the B^j, the C^j, Q, and b
    BB <- lapply( AA, "%*%", hts )
    CC <- sapply( BB, function (B) { B[all.zeros] <- 0; rowSums(B) } )
    b <- (-1) * colSums(CC) * ncol(hts)
    Q <- crossprod(CC)
    # estimates of aa!!
    return( solve( Q, b ) )
}

check.aa <- function (aa,hts) {
    # check for consistency
    G <- aa[1] * AA[[1]]
    for (k in 2:n.layers) {
        G <- G + aa[k] * AA[[k]]
    }
    ones <- (-1) * G %*% hts
    range(ones[-all.zeros])
}

# check this works with truth
estimate.aa(true.hts)
# consistent?
check.aa(.Last.value,true.hts)


# and interpolation
est.aa <- estimate.aa(interp.hts)
check.aa(est.aa,true.hts)
check.aa(est.aa,interp.hts)
