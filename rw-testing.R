#!/usr/bin/Rscript

# Testing code

source("rw-testing-fns.R")

####
# Random walk on a directed graph:
#   walk moves at rate 1
#   to a randomly chosen location
#   out of the possibilities

nlocs <- 4
adj <- matrix(rbinom(nlocs^2,size=1,prob=3/nlocs), nrow=nlocs, ncol=nlocs)
while ( any( rowSums( adj ) == 0 ) || any( colSums( adj ) == 0 ) ) { 
    adj <- matrix(rbinom(nlocs^2,size=1,prob=3/nlocs), nrow=nlocs, ncol=nlocs)
}
diag(adj) <- 0
adj <- ( adj > 0 )

nwalks <- 5000
nsteps <- 200
xt <- simrw(x0=sample.int(nlocs,size=nwalks,replace=TRUE),t=nsteps,adj=adj)

# hitting times
htimes <- hitting(xt,locs=1:nlocs)
# matrix of mean hitting times: hmat[x,y] is hitting time of y from x
hmeans <- do.call( rbind, tapply( 1:ncol(htimes), xt[1,], function (k) { rowMeans(htimes[,k,drop=FALSE],na.rm=TRUE) } ) )

# generator matrix
gen <- adj.to.gen(adj)

# if ty[x] is mean hitting time of y beginning at x, with ty[y]=0, then should have:
#   (-gen[-y,-y]) %*% ty[-y] = 1
onesq <- sapply( 1:nlocs, function (k) { (-1) * gen[-k,-k] %*% hmeans[-k,k] } )

# these are close to 1!
onesq
