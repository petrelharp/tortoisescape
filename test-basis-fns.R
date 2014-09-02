#!/usr/bin/R

# see
#  http://lmb.informatik.uni-freiburg.de/papers/download/wa_report01_08.pdf?origin=publication_detail
#  ( http://www.bibsonomy.org/bibtex/2e3ede022c042c956308ca4055b7af1ef/peter.ralph )


source("rw-testing-fns.R")
source("basis-fns.R")

# width of square grid
n <- 1001
nmap <- matrix(1:n^2,nrow=n)
all.locs <- data.frame( i=as.vector(row(nmap)), j=as.vector(col(nmap)) )

# random landscape
G <- grid.generator(n)
orig <- ij.to.k(floor(c(n/2,n/2)),n)
Gz <- G[-orig,-orig]


##
# analytic:

# check 'em out: Bessel
xx <- seq(0,10,length.out=400)
nfuns <- 100
B.zeros <- bessel_zero_J1(1:nfuns)
B.a <- 10
Bfuns <- bessel_J1( outer( xx, B.zeros/B.a, "*" ) )
dim(Bfuns) <- c(length(xx),length(B.zeros))
matplot( matrix(xx[row(Bfuns)],nrow=nrow(Bfuns)), Bfuns, type='l')


# bivariate
xy <- as.matrix(all.locs)/n

k <- 2; m <- 5
radfun <- radial(xy,k=k,m=m,x0=xy[orig,])
argfun <- angular(xy,m=m,x0=xy[orig,])
dim(radfun) <- dim(argfun) <- c(n,n)

image(Matrix(radfun))
image(Matrix(argfun))
image(Matrix(argfun*radfun))


##
# Empirical approach:

# classical:
if (n<102) {
    eGz <- eigen(Gz)

    for (k in 1:100) {
        tmp <- eGz$vectors[,k]
        tmp <- c(tmp[1:(orig-1)],0,tmp[orig:length(tmp)])
        dim(tmp) <- c(n,n)
        plot( image(Matrix(tmp)) )
        readline("next?")
    }

}

# ARPACK
# and https://github.com/yixuan/rARPACK
require(rARPACK)

layout(matrix(1:12,nrow=3))
for (k in 1:12) {
    # with( all.locs, plot( i, j, col=colorize(v1$vectors[,k]) ), pch=20, cex=.25 )
    bigones <- which( abs(v1$vectors[,k]) > .01 )
    with( all.locs, plot( i[bigones], j[bigones], col=heat.colors(32)[cut(v1$vectors[,k],breaks=seq(-1,1,length.out=24))], pch=20, xlim=range(i), ylim=range(j) ) )
}


# or: see http://cran.r-project.org/web/packages/svd/index.html
require(svd)

Gprod <- function (x) { as.vector( Gz %*% x ) }
tGprod <- function (x) { as.vector(crossprod( Gz, x )) }
ext.Gz <- extmat( Gprod, tGprod, nrow=nrow(Gz), ncol=ncol(Gz) )

system.time( v1 <- propack.svd(ext.Gz, neig = 10) )
system.time( v2 <- trlan.svd(ext.Gz, neig = 10) )

plot( v1$u[,1], Gprod(v1$u[,1]) )
abline(0,v1$d[1])

