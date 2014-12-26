#!/usr/bin/Rscript

source("resistance-fns.R")

# width of square grid
n <- 100
nmap <- matrix(1:n^2,nrow=n)
all.locs <- cbind( i=as.vector(row(nmap)), j=as.vector(col(nmap)) )

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


# Observed hitting times
pairwise.hts <- true.hts[locs,]



###########
# now try to interpolate just one of these
#    quadratic program

# with varying noise values

kk <- 3
obs.hts <- pairwise.hts[,kk]
epsvals <- c(0,10^(seq(-4,-1,length.out=10)))
interp.hts.list <- lapply( epsvals, function (eps) {
        noisy.obs.hts <- obs.hts * exp( eps * rnorm(length(obs.hts)) )
        gamma <- 1
        Pmat <- sparseMatrix( i=seq_along(locs), j=locs, x=1, dims=c(length(locs),nrow(G)) )
        PtP <- gamma * crossprod(Pmat)
        GtG <- crossprod( G[-locs[kk],] )
        bvec <- gamma * crossprod(Pmat,noisy.obs.hts) - crossprod( G[-locs[kk],], rep(1.0,nrow(G)-1) )
        # bvec <- gamma * crossprod(Pmat,noisy.obs.hts) + rexp(length(bvec))
        # bvec <- rexp(length(bvec))
        interp.hts <- solve( PtP+GtG, bvec )
        return(interp.hts)
    } )

layout(matrix(c(1,2,1,3),nrow=2))
for (k in seq_along(interp.hts.list)) {
    plot( interp.hts.list[[k]], true.hts[,kk], main=epsvals[k] ); abline(0,1)
    tmp <- true.hts[,kk]; dim(tmp) <- c(n,n); image(tmp)
    tmp <- as.numeric(interp.hts.list[[k]]); dim(tmp) <- c(n,n); image(tmp)
    if (is.null(locator(1))) { break }
}

###
# all of them
all.interp.hts <- interp.hitting( G, locs, pairwise.hts, locs )

range( all.interp.hts - true.hts )

for (kk in seq_along(locs)) {
    layout(matrix(c(1,2,1,3),nrow=2))
    plot( all.interp.hts[,kk], true.hts[,kk] )
    tmp <- true.hts[,kk]; dim(tmp) <- c(n,n); image(tmp)
    tmp <- as.numeric(all.interp.hts[,kk]); dim(tmp) <- c(n,n); image(tmp)
    if (is.null(locator(1))) { break }
}
