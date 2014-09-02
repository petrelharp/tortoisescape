#!/usr/bin/Rscript

require(Matrix)
source("resistance-fns.R")

n <- 100

b <- rexp(n^2)
A <- grid.generator(n,killing=rexp(n^2))

# n=1000 -> 0.016s
system.time( Ab <- A %*% b )


# n=1000 -> 510.780s for solve but only 0.016s for %*%
system.time( x.exact <- solve(A,b) )


# try conjugate gradient?
require(Rcgmin)

AtA <- crossprod(A)
Atb <- crossprod(A,b)
f <- function (x) { sum( ( A%*%x - b )^2 ) }
g <- function (x) { 2 * ( AtA %*% x - Atb ) }

x0 <- b

# test gradient
f0 <- f(x0)
kk <- sample.int(length(x0),100)
test.gr <- sapply( kk, function (k) { x1 <- x0; x1[k] <- x1[k]+.001; 1000 * ( f(x1) - f0 ) } )
range( (test.gr - g(x0)[kk])/abs(g(x0)[kk]) )  # yup

system.time( est.x <- Rcgmin( x0, fn=f, gr=g, control=list(trace=1,maxit=2000) ) )
# system.time( est.x.1 <- Rcgmin( est.x$par, fn=f, gr=g, control=list(trace=1) ) )
