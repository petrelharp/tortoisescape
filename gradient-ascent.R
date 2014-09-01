#!/usr/bin/Rscript

source("rw-testing-fns.R")

# width of square grid
n <- 30
nmap <- matrix(1:n^2,nrow=n)
all.locs <- data.frame( i=as.vector(row(nmap)), j=as.vector(col(nmap)) )

# random landscape
n.layers <- 2
AA <- replicate( n.layers, grid.generator(n), simplify=FALSE )
aa <- rep(1/n.layers,length=n.layers)
make.G <- function (aa) {
    G <- aa[1] * AA[[1]]
    for (k in 2:n.layers) {
        G <- G + aa[k] * AA[[k]]
    }
    return(G)
}
G <- make.G(aa)

# sampling locations
nsamps <- 20
locs.ij <- data.frame( i=sample.int(n,nsamps)-1L, j=sample.int(n,nsamps)-1L )
locs <- ij.to.k(locs.ij,n)

# analytical mean hitting times
true.hts <- hitting.analytic(locs,G)  # warning, this could take a while (10s for n=100 and nsamps=20)

# Observed hitting times
pairwise.hts <- true.hts[locs,]

L <- function (aa) {
    # objective function: alphas
    assign("G", make.G(aa), parent.env(environment()))
    assign("hts", hitting.analytic(locs,G), parent.env(environment()))
    sum( ( hts[locs,] - pairwise.hts )^2 )
}
dL <- function (aa) {
    # gradient of L, wrt aa:
    sapply( seq_along(aa), function (i) {
            AH <- AA[[i]] %*% hts
            GAH <- sapply( seq_along(locs), function (k) as.vector( solve( G[-locs[k],-locs[k]], AH[-locs[k],k] ) ) )
            return( -2 * sum( GAH[locs,] * (hts[locs,]-pairwise.hts) ) )
        } )
}
environment(L) <- environment(dL) <- fun.env <- list2env( list(G=G,hts=true.hts), parent=environment() )


# WORKS: with 250 function evaluations
optim( par=aa+.1*runif(2,min=-1,max=1), fn=L, gr=dL, control=list(parscale=c(.1,.1),fnscale=1e5,trace=5), method="CG" )

# WORKS: with 500 function evaluations
optim( par=c(.2,1), fn=L, gr=dL, control=list(parscale=c(.1,.1),fnscale=1e5,trace=5), method="CG" )


if (FALSE) {

    ## CHECKING THE MATH:

    # check the gradient
    a0 <- c(.49,.54)
    eps <- c(.00001,-.00001)
    tmp <- c( L(a0+eps), L(a0-eps), L(a0) )
    c( sum(dL(a0)*eps), (tmp[1]-tmp[3]), (tmp[3]-tmp[2]) ) # YEP, well, relative error like 1e-3

    # check derivative of H:
    i <- 1; k <- 1
    eps <- .001
    a0 <- c(.48,.51)
    gg <- make.G(c(.48,.51))
    ah <- AA[[i]] %*% hitting.analytic( locs[k], gg )
    gah <- as.vector( solve( gg[-locs[k],-locs[k]], ah[-locs[k],1] ) )
    rgah <- ( hitting.analytic( locs[k], gg ) - hitting.analytic( locs[k], make.G(a0+c(eps,0)) ) )[-locs[k]]
    summary( (rgah - eps * gah)/gah ) # YUP of order 1e-6


####
# the following seems faster but is not

pairwise.hts.nozeros <- sapply( seq_along(locs), function (k) { pairwise.hts[-k,k] } )
locs.nozeros <- cbind( locs, rep(seq_along(locs),each=length(locs)) )
locs.nozeros <- locs.nozeros[ match(locs.nozeros[,1],locs) != locs.nozeros[,2], ]
L1 <- function (aa) {
    # objective function: alphas
    G <- make.G(aa)
    hts <- sapply( seq_along(locs), function (k) {
            solve( G[-locs[k],-locs[k]], rep.int(-1.0,nrow(G)-1L) )
        } )
    sum( ( hts[locs.nozeros] - pairwise.hts.nozeros )^2 )
}
dL1 <- function (aa) {
    # gradient of L, wrt aa:
    sapply( seq_along(aa), function (i) {
            AH <- AA[[i]] %*% hts
            GAH <- sapply( seq_along(locs), function (k) solve( G[-locs[k],-locs[k]], AH[-k] ) )
            return( rowSums( GAH[locs.nozeros] * (hts[locs.nozeros]-pairwise.hts.nozeros) ) )
        } )
}
environment(L1) <- environment(dL1) <- fun.env1 <- list2env( list(G=G,hts=true.hts), parent=environment() )

system.time( L1(aa) )
system.time( dL1(aa) )

#########
# testing:  can pass something back and forth between function and gradient in optim?
#  yes!

x <- matrix( rnorm(200), nrow=2 )
sig <- 1

f <- function (x0) {
    assign( "sig", sqrt( mean( (x-x0)^2 ) ), env=fenv )
    return( sum( (x-x0)^2/sig^2 ) )
}
fenv <- environment(f)

g <- function (x0) {
    return( 2 * rowSums( x-x0 ) / sig^2 )
}
environment(g) <- fenv

optim( par=c(0.1,0.1), fn=f, gr=g, control=list(fnscale=-1), method="CG" )

g(0) == 2 * rowSums( x )  # FALSE!

}
