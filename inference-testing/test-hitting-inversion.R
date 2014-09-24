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

analytic.true.hts <- hitting.analytic(locs,G)  # warning, this could take a while (10s for n=100 and nsamps=20)

# jacobi wants a G WITHOUT a diagonal
jG <- G
diag(jG) <- 0

hts <- matrix(rnorm(nrow(G)*length(locs)),nrow=nrow(G),ncol=length(locs))

jacobi.true.hts <- hitting.jacobi(locs,jG,analytic.true.hts)
jacobi.true.hts <- hitting.jacobi(locs,jG,ifelse(hts>0,1,0))

layout(1:2)
for (k in 1:ncol(hts)) {
    image( matrix(analytic.true.hts[,k],nrow=length(locs)) )
    image( matrix(jacobi.true.hts[,k],nrow=length(locs)) )
    if (is.null(locator(1))) { break }
}

plot(analytic.true.hts[,1],jacobi.true.hts[,1])


