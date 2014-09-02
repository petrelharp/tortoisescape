#!/usr/bin/Rscript

source("rw-testing-fns.R")
source("basis-fns.R")

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
locs <- ij.to.k(locs.ij,n)

# analytical mean hitting times
true.hts <- hitting.analytic(locs,G)  # warning, this could take a while (10s for n=100 and nsamps=20)


# Observed hitting times
pairwise.hts <- true.hts[locs,]

# now try to interpolate just one of these
kk <- 1
obs.hts <- pairwise.hts[,kk]


# first-order radial smooth
all.rads <- rowSums( sweep((all.locs),2,locs.ij[kk,],"-")^2 )
rads <- all.rads[locs]
trads <- rads^(1/4)
ht.smooth <- loess( obs.hts ~ trads, span=.5, control=loess.control(surface='direct') )
hts.interp.0 <- pmax(0, predict( ht.smooth, newdata=data.frame(trads=all.rads^(1/4)) ) )

layout(matrix(c(1,2,1,3),nrow=2))
plot( true.hts[,kk], hts.interp.0, col=ifelse(1:nrow(true.hts) %in% locs, 'red','black'), pch=20 )
abline(0,1)
tmp <- true.hts[,kk]; dim(tmp) <- c(n,n); image(tmp)
tmp <- hts.interp.0; dim(tmp) <- c(n,n); image(tmp)


# build the basis
xy <- as.matrix(all.locs)/n
fbasis <- fbbasis(xy,mmax=30,x0=xy[locs[kk],])
fbasis.gram <- crossprod(fbasis)

# check
for (ii in 1:ncol(fbasis)) {
    tmp <- matrix( fbasis[,ii],nrow=n,ncol=n )
    plot( image(Matrix(tmp)) )
    # plot( j ~ i, data=all.locs, bg=adjustcolor(ifelse(fbasis[,ii]>0,"blue","red"),.25), pch=21, cex=scale(abs(fbasis[,ii]))*3 )
    readline("next?")
}

# orthogonal-ish?
image(Matrix(crossprod(cbind( hts.interp.0/sqrt(sum(hts.interp.0^2)),fbasis))))

# project into these coords
require(MASS)
resid.hts <- (obs.hts - hts.interp.0)
resid.hts.proj <- rowSums( sweep( fbasis, 2, ginv(fbasis.gram) %*% crossprod(fbasis,resid.hts), "*" ) )

hts.interp <- hts.interp.0 + resid.hts.proj

plot( true.hts[,kk], hts.interp )

layout(t(1:2))
tmp <- hts.interp; dim(tmp) <- c(n,n); image(tmp)
tmp <- true.hts[,kk]; dim(tmp) <- c(n,n); image(tmp)
