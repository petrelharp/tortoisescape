#!/usr/bin/Rscript


# mappings between index in a matrix of height (in rows) n
#   ZERO-BASED: (i,j) and column-oriented (k)
#   ALLOW indices outside the grid
.ob <- function (ij,n,m){ ( ij[,1] >= 0 ) & ( ij[,1] < n ) & ( ij[,2] >= 0 ) & ( ij[,2] < m ) }
ij.to.k <- function (ij,n,m) { if (is.null(dim(ij))) { dim(ij) <- c(1,length(ij)) }; ifelse( .ob(ij,n,m), ij[,1,drop=FALSE]+ij[,2,drop=FALSE]*n, NA ) }
k.to.ij <- function (k,n) { cbind( k%%n, k%/%n ) }
shift <- function (dij,k,n,m) { ij.to.k( sweep( k.to.ij(k,n), 2, as.integer(dij), "+" ), n, m ) }

require(Matrix)

grid.dist <- function (ij) {
    # return matrix of pairwise as-the-crow-flies distances
    sqrt( outer(ij[,1],ij[,1],"-")^2 + outer(ij[,2],ij[,2],"-")^2 )
}

grid.adjacency <- function (n,m=n,diag=TRUE,symmetric=TRUE) {
    # for a grid of height n and width m
    nn <- 0:(n*m-1)
    nreps <- if(diag){5}else{4}
    adj <- data.frame(
            i=rep(nn, nreps),
            j=c( if(diag){nn}else{NULL},
                 shift(c(+1,0),nn,n,m),
                 shift(c(-1,0),nn,n,m),
                 shift(c(0,+1),nn,n,m),
                 shift(c(0,-1),nn,n,m)
                 ),
            x=rep(1,nreps*length(nn))
            )
    # on the boundary?
    usethese <- ! ( is.na(adj$j) | !.ob(k.to.ij(adj$j,n),n,m) | !.ob(k.to.ij(adj$i,n),n,m) )
    if (symmetric) { usethese <- ( usethese & ( adj$i < adj$j ) ) }
    # add 1 since we worked 0-based above; sparseMatrix (unlike underlying representation) is 1-based.
    A <- with( subset(adj, usethese ), sparseMatrix( i=i+1L, j=j+1L, x=x, dims=c(n*m,n*m), symmetric=symmetric ) )
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
                z <- solve( G[-k,-k], rep.int(-1.0,nrow(G)-1L) ) 
                return( c( z[seq.int(1,length.out=k-1L)], 0, z[seq.int(k,length.out=length(z)-k+1L)] ) )
            } )
    return(hts)
}

interp.hitting <- function ( G, locs, obs.hts, gamma=1 ) {
    # interpolate hitting times by minimizing squared error:
    #       G is a generator matrix
    #       locs is the indices of the rows of G for which we have data
    #       obs.hts is the (locs x locs) matrix of mean hitting times
    #       gamma is a fudge factor (so far unnecessary?)
    #   note that without the second 'crossprod' term in the defn of 'bvec'
    #   this finds the harmonic function interpolating the observed hitting times
    #   (or something close, I think)
    Pmat <- sparseMatrix( i=seq_along(locs), j=locs, x=1, dims=c(length(locs),nrow(G)) )
    PtP <- gamma * crossprod(Pmat)
    sapply( seq_along(locs), function (kk) {
                Gk <- G[-locs[kk],]
                bvec <- gamma * crossprod(Pmat,obs.hts[,kk]) - crossprod( Gk, rep(1.0,nrow(G)-1) )
                as.numeric( solve( PtP+crossprod(Gk), bvec ) )
    } )
}

make.G <- function (aa,AA) {
    G <- aa[1] * AA[[1]]
    if (length(AA)>1) for (k in 2:length(AA)) {
        G <- G + aa[k] * AA[[k]]
    }
    return(G)
}

estimate.aa <- function (hts,locs,AA) {
    # Estimate alphas given full hitting times
    # don't count these cells:
    zeros <- locs + (0:(length(locs)-1))*nrow(hts)
    # here are the B^j, the C^j, Q, and b
    BB <- lapply( AA, "%*%", hts )
    CC <- sapply( BB, function (B) { B[zeros] <- 0; rowSums(B) } )
    b <- (-1) * colSums(CC) * ncol(hts)  # WHY THIS EXTRA FACTOR OF m?
    Q <- crossprod(CC)
    return( solve( Q, b ) )
}

iterate.aa <- function (aa,hts,locs,AA) {
    G <- make.G(aa,AA)
    # interpolate hitting times
    interp.hts <- interp.hitting( G, locs, hts )
    # infer aa from these
    estimate.aa(interp.hts,locs,AA)
}

##
# plotting whatnot

plot.ht <- function (ht,dims=c(sqrt(length(ht)),sqrt(length(ht)))) {
    dim(ht) <- dims; image(ht)
}

plot.hts <- function (hts,dims=c(sqrt(nrow(hts)),sqrt(nrow(hts)))) {
    for (k in 1:ncol(hts)) {
        plot.ht(hts[,k],dims=dims)
        readline("next?")
    }
}

colorize <- function (x, nc=32, colfn=function (n) rainbow_hcl(n,c=100,l=50), zero=FALSE, trim=0, breaks, return.breaks=FALSE) {
    if (is.numeric(x) & trim>0) {
        x[ x<quantile(x,trim,na.rm=TRUE) ] <- quantile(x,trim,na.rm=TRUE)
        x[ x>quantile(x,1-trim,na.rm=TRUE) ] <- quantile(x,1-trim,na.rm=TRUE)
    }
    if (missing(breaks) & is.numeric(x)) {
        if (zero) {
            breaks <- seq( (-1)*max(abs(x),na.rm=TRUE), max(abs(x),na.rm=TRUE), length.out=nc )
        } else {
            breaks <- seq( min(x,na.rm=TRUE), max(x,na.rm=TRUE), length.out=nc )
        }
        x <- cut(x,breaks=breaks,include.lowest=TRUE)
    } else {
        x <- factor(x)
    }
    if (return.breaks) {
        return(breaks)
    } else {
        return( colfn(nlevels(x))[as.numeric(x)] )
    }
}
