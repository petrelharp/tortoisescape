#!/usr/bin/Rscript

# http://stackoverflow.com/questions/7547758/using-sample-with-sample-space-size-1
resample <- function(x, ...) x[sample.int(length(x), ...)]

simrw <- function (x0,t,adj) {
    # sample RW on graph defined by adj
    xx <- matrix(NA,nrow=t+1,ncol=length(x0))
    xx[1,] <- x0
    for (s in 1:t) {
        steps <- lapply( xx[s,], function (k) which(adj[k,]) )
        xx[s+1,] <- sapply( steps, resample, size=1 )
    }
    return(xx)
}

plot.rw <- function (rw,locs,...) {
    visited <- unique(unlist(rw))
    t <- nrow(rw)-1
    plot( locs[visited,], xlim=range(locs$x), ylim=range(locs$y), ... )
    for (k in 1:ncol(rw)) {
        segments( x0=locs$x[rw[-(t+1),k]], y0=locs$y[rw[-(t+1),k]], x1=locs$x[rw[-1,k]], y1=locs$y[rw[-1,k]], col=rainbow(ncol(rw)+4)[k], lwd=2 )
    }
}

coalrw <- function (rw,coalprob=1) {
    # sample coalescent times given paths of random walks
    # IGNORE dependencies for now -- does NOT get JOINT coalescent times
    ctimes <- data.frame( t( combn( 1:ncol(rw), 2 ) ) )
    ctimes$ctime <- NA
    csteps <- 1+rgeom( 1:nrow(ctimes), coalprob )  # coalesces after this many steps together
    for (k in 1:nrow(ctimes)) {
        cotimes <- which( rw[,ctimes[k,1]]==rw[,ctimes[k,2]] )
        if ( length(cotimes) >= csteps[k] ) {
            ctimes$ctime[k] <- cotimes[csteps[k]]
        }
    }
    return(ctimes)
}

hitting <- function (rw,locs=sort(unique(as.vector(rw))) ) {
    # for each sample path in rw
    # find the vector of hitting times of each location
    apply( rw, 2, function (x) { sapply( locs, function (loc) { match(loc,x)-1 } ) } )
}

adj.to.gen <- function (adj) {
    # From adjacency matrix return generator for:
    # random walk on a directed graph that
    #   walk moves at rate 1
    #   to a randomly chosen location
    #   out of the possibilities
    gen <- sweep(adj,1,rowSums(adj),"/")
    diag(gen) <- (-1)*rowSums(gen)
    return(gen)
}


# mappings between index in a square matrix
#   ZERO-BASED: (i,j) and column-oriented (k)
#   ALLOW indices outside the grid
#   grid height (number of rows) is n
.ob <- function (ij,n){ ( ij[,1] >= 0 ) & ( ij[,1] < n ) & ( ij[,2] >= 0 ) & ( ij[,2] < n ) }
ij.to.k <- function (ij,n) { ifelse( .ob(ij,n), ij[,1,drop=FALSE]+ij[,2,drop=FALSE]*n, NA ) }
k.to.ij <- function (k,n) { cbind( k%%n, k%/%n ) }
shift <- function (dij,k,n) { ij.to.k( sweep( k.to.ij(k,n), 2, as.integer(dij), "+" ), n ) }

require(Matrix)

grid.adjacency <- function (n) {
    nn <- 0:(n^2-1)
    adj <- data.frame(
            i=rep(nn,5),
            j=c( nn,
                 shift(c(+1,0),nn,n),
                 shift(c(-1,0),nn,n),
                 shift(c(0,+1),nn,n),
                 shift(c(0,-1),nn,n)
                 ),
            x=rep(1,5*length(nn))
            )
    boundary <- is.na(adj$j) | !.ob(k.to.ij(adj$j,n),n) | !.ob(k.to.ij(adj$i,n),n)
    utri <- ( adj$i < adj$j )
    A <- with( subset(adj,!boundary & utri ), sparseMatrix( i=i+1L, j=j+1L, x=x, dims=c(n^2,n^2), symmetric=TRUE ) )
    # A <- as( with( subset(adj,!boundary), new( "dgTMatrix", i=i, j=j, x=x, Dim=as.integer(c(n^2,n^2)) ) ), "dgCMatrix" )
    return(A)
}

grid.generator <- function (n,killing=0) {
    # random generator for RW with killing on square 2D grid
    A <- grid.adjacency(n)
    A@x <- rexp(length(A@x))
    A <- ( A + t(A) )
    diag(A) <- (-1)*rowSums( A - Diagonal(nrow(A),diag(A)) ) - killing
    return(A)
}

hitting.analytic <- function (locs,G) {
    # compute analytical expected hitting times
    hts <- sapply( locs, function (k) { 
                z <- solve( G[-k,-k], rep(-1,nrow(G)-1L) ) 
                return( c( z[1:(k-1)], 0, z[k:length(z)] ) )
            } )
    return(hts)
}

##
# plotting whatnot

colorize <- function (x, nc=32, colfn=function (n) rainbow_hcl(n,c=100,l=50), zero=FALSE, trim=0) {
    require(colorspace)
    if (is.numeric(x) & trim>0) {
        x[ x<quantile(x,trim,na.rm=TRUE) ] <- quantile(x,trim,na.rm=TRUE)
        x[ x>quantile(x,1-trim,na.rm=TRUE) ] <- quantile(x,1-trim,na.rm=TRUE)
    }
    if (is.numeric(x)) {
        if (zero) {
            breaks <- seq( (-1)*max(abs(x),na.rm=TRUE), max(abs(x),na.rm=TRUE), length.out=nc )
        } else {
            breaks <- seq( min(x,na.rm=TRUE), max(x,na.rm=TRUE), length.out=nc )
        }
        x <- cut(x,breaks=breaks,include.lowest=TRUE)
    } else {
        x <- factor(x)
    }
    return( colfn(nlevels(x))[as.numeric(x)] )
}

