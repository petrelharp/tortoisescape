#!/usr/bin/R

####
# Comparing various matrix calculations
#
# Conclusion: coercing vectors to Matrix helps, a bit.

require(Matrix)
require(microbenchmark)

source("resistance-fns.R")

n <- 100

A <- grid.generator(n,killing=rexp(n^2))
AC <- as(A,"dgCMatrix")   # doesn't hurt

b <- rexp(n^2)
bM <- Matrix(b)
bC <- as(bM,"dgCMatrix")

microbenchmark( A%*%b )
microbenchmark( A%*%bM )  # best
microbenchmark( A%*%bC )

microbenchmark( AC%*%b )
microbenchmark( AC%*%bM )
microbenchmark( AC%*%bC )
