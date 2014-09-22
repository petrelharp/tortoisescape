#!/usr/bin/Rscript

nn <- floor( exp( log(10)/4 * 1:(4*3.5) ) )
nreps <- max(nn)/nn

A <- matrix( rnorm( max(nn)^2 ), nrow=max(nn) )
b <- rnorm( max(nn) )

nn.times <- sapply( seq_along(nn), function (k) { n <- nn[k]; An <- A[1:n,1:n]; bn <- b[1:n]; system.time( for (k in 1:nreps[k]) solve(An,bn) ) } )

pdf(file="dense-inversion-timing.pdf")
plot( nn, nn.times["user.self",]/nreps, log='xy', xlab="dimension of matrix", ylab="seconds to solve( )" )
dev.off()
