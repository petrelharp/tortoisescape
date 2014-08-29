#!/usr/bin/Rscript

require(Matrix)
source("rw-testing-fns.R")

# results
# n=1000 -> 510.780s for solve but only 0.016s for %*%

n <- 100

b <- rexp(n^2)
A <- grid.generator(n,killing=rexp(n^2))

system.time( Ab <- A %*% b )

system.time( x <- solve(A,b) )

luA <- lu(A)  # these factors are 33x less sparse than A (but still sparse...)
image(luA@L)
dev.set(dev.next())
image(luA@U)
