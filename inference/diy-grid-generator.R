
# mappings between index in a matrix of height (in rows) n
#   ZERO-BASED: (i,j) and column-oriented (k)
#   ALLOW indices outside the grid
.ob <- function (ij,n,m){ ( ij[,1] >= 0 ) & ( ij[,1] < n ) & ( ij[,2] >= 0 ) & ( ij[,2] < m ) } # outside grid
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
                 )
            )
    # on the boundary?
    usethese <- ! ( is.na(adj$j) | !.ob(k.to.ij(adj$j,n),n,m) | !.ob(k.to.ij(adj$i,n),n,m) )
    if (symmetric) { usethese <- ( usethese & ( adj$i <= adj$j ) ) }
    # add 1 since we worked 0-based above; sparseMatrix (unlike underlying representation) is 1-based.
    A <- with( subset(adj, usethese ), sparseMatrix( i=i+1L, j=j+1L, x=1.0, dims=c(n*m,n*m), symmetric=symmetric ) )
    return(A)
}

grid.sum <- function (n,m=n,direction,zeros=FALSE) {
    # for a grid of height n and width m,
    #   return operator that sums each site with the other one at displacement 'direction'
    # if there is no such site, either:
    # if (zeros) : include zero rows
    # else : omit these,
    #   and SHOULD RETURN IN SAME ORDER AS grid.adjacency
    nn <- 0:(n*m-1)
    nreps <- 2
    adj <- data.frame(
            i=rep(nn, nreps),
            j=c( nn,
                 shift(direction,nn,n,m)
                 ),
            x=rep(1,nreps*length(nn))
            )
    # on the boundary?
    usethese <- with( adj[n*m+(1:(n*m)),], ! ( is.na(j) | !.ob(k.to.ij(j,n),n,m) | !.ob(k.to.ij(i,n),n,m) ) )
    adj <- subset( adj, c(usethese,usethese) )
    if (!zeros) {
        # remove zero-rows
        zrows <- cumsum( tabulate( 1L+adj$i, nbins=n*m ) > 0 )
        adj$i <- zrows[ 1L+adj$i ] - 1L
    }
    # add 1 since we worked 0-based above; sparseMatrix (unlike underlying representation) is 1-based.
    A <- with( adj, sparseMatrix( i=i+1L, j=j+1L, x=x, dims=c(if(zeros){n*m}else{zrows[n*m]},n*m) ) )
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

###
# OLD STUFF
#   ...unused except in older scripts


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

