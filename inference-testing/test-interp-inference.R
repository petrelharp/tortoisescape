#!/usr/bin/Rscript
##
# Build the map of which sets of alpha parameters map to each other

source("resistance-fns.R")
require(parallel)

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

###
# get grid of where parameters map to:
# TAKES A LONG TIME
avals <- (seq(.01,sqrt(5),length.out=40))^2
agrid <- as.matrix( expand.grid( a0=avals, a1=avals ) )
bgrid <- mclapply( 1:nrow(agrid), function (k) {
        a0 <- agrid[k,]
        iterate.aa(a0,pairwise.hts,locs,AA) }, mc.cores=16 )
bgrid <- do.call( rbind, bgrid )
colnames(bgrid) <- c("b0","b1")

pdf(file="abmap.pdf",pointsize=10,width=10,height=10)
plot( 0, 0, xlim=range(agrid[,1]), ylim=range(agrid[,2]), type='n', xlab='', ylab='' )
arrows( x0=agrid[,1], y0=agrid[,2], x1=bgrid[,1], y1=bgrid[,2], length=0.1 )
points( aa[1], aa[2], pch="*", col="red", cex=2 )
dev.off()

write.table( cbind(agrid,bgrid), file="abgrid.csv", sep=",", row.names=FALSE )


