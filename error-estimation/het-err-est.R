

###
# binomial error
nsites <- 1e6

phet <- .1
perr <- .01
mean.coverage <- 1

coverage <- rpois(nsites,mean.coverage)
het <- ( rbinom(nsites,size=1,prob=phet) == 1 )
probs <- ifelse( het, 1/2, perr )
counts <- rbinom( nsites, size=coverage, prob=probs )
counts.ref <- coverage-counts

likfun <- function (pp) {
    # pp is ( phet, perr )
    sum( log( pp[1] * 2^(-coverage) + (1-pp[1]) * pp[2]^counts * (1-pp[2])^counts.ref ) )
}

pvals <- expand.grid( phet=seq(.8*phet,1.2*phet,length.out=20), perr=seq(.8*perr,1.2*perr,length.out=20) )
pvals$lik <- apply( pvals, 1, likfun )

cplot <- function (x,y,z,...) {
    xvals <- sort(unique(x))
    yvals <- sort(unique(y))
    zmat <- matrix( z, nrow=length(xvals), ncol=length(yvals) )
    image(xvals,yvals,zmat,...)
    contour(xvals,yvals,zmat,...,add=TRUE)
}

with( pvals, cplot( phet, perr, lik ) )
abline(v=phet,h=perr)

layout(t(1:2))
with(pvals, plot( phet, lik, col=cut(perr,8) ) )
abline(v=phet)
with( subset(pvals,abs(perr-get("perr",parent.env(environment())))==min(abs(perr-get("perr",parent.env(environment()))))), lines(phet, lik ) )
with(pvals, plot( perr, lik, col=cut(phet,8) ) )
abline(v=perr)



###
# beta-distributed error
nsites <- 1e6

phet <- .1
perr <- .01
nerr <- 10
mean.coverage <- 10

coverage <- rpois(nsites,mean.coverage)
het <- ( rbinom(nsites,size=1,prob=phet) == 1 )
probs <- ifelse( het, 1/2, rbeta( nsites, shape1=perr*nerr, shape2=(1-perr)*nerr ) )
counts <- rbinom( nsites, size=coverage, prob=probs )
counts.ref <- coverage-counts

likfun <- function (pp) {
    # pp is ( phet, perr, nerr )
    sum( log( pp[1] * 2^(-coverage) + (1-pp[1]) * beta(counts+pp[2]*pp[3],counts.ref+(1-pp[2])*pp[3])/beta(pp[2]*pp[3],(1-pp[2]*pp[3])) ) )
}

pvals <- expand.grid( phet=seq(.8*phet,1.2*phet,length.out=20), perr=seq(.8*perr,1.2*perr,length.out=20), nerr=10 )
pvals$lik <- apply( pvals, 1, likfun )

cplot <- function (x,y,z,...) {
    xvals <- sort(unique(x))
    yvals <- sort(unique(y))
    zmat <- matrix( z, nrow=length(xvals), ncol=length(yvals) )
    image(xvals,yvals,zmat,...)
    contour(xvals,yvals,zmat,...,add=TRUE)
}

with( pvals, cplot( phet, perr, lik ) )

layout(t(1:2))
with(pvals, plot( phet, lik, col=cut(perr,8) ) )
with(pvals, plot( perr, lik, col=cut(phet,8) ) )
