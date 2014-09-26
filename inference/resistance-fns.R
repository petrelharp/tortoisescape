#!/usr/bin/Rscript

# move @p to @j in dgCMatrix
p.to.j <- function (p) { rep( seq.int( length(p)-1 ), diff(p) ) }

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

hitting.jacobi <- function (locs,G,hts,idG=1/rowSums(G),b=-1.0,tol=1e-6,kmax=1000) {
    # compute analytical expected hitting times using the Jacobi method (from jacobi.R)
    #  note that G ** comes with no diagonal **
    for (locnum in 1:length(locs)) {
        ll <- locs[locnum]
        k <- 1
        x <- hts[,locnum]
        x[ll] <- 0
        for (k in 1:kmax) {
            x_new <- idG * (G%*%x-b)
            x_new[ll] <- 0
            err <- mean((x_new-x)^2)
            # cat(k,":",err,"\n")
            if (err < tol) {
                cat("converged! err=", err, "\n")
                break; 
            }
            x <- x_new
        }
        if (k==kmax) { cat("Hit kmax. Did not converge, err=",err,"\n") }
        hts[,locnum] <- as.vector(x_new)
    }
    return(hts)
}

hitting.analytic <- function (locs,G) {
    # compute analytical expected hitting times
    if ("parallel" %in% .packages()) {
        numcores<-as.numeric(scan(pipe("cat /proc/cpuinfo | grep processor | tail -n 1 | awk '{print $3}'")))+1
        this.apply <- function (...) { do.call( cbind, mclapply( ..., mc.cores=numcores ) ) }
    } else {
        this.apply <- function (...) { sapply( ... ) }
    }
    hts <- this.apply( locs, function (k) { 
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
# misc

selfname <- function (x) { names(x) <- make.names(x); x }

##
# plotting whatnot

plot.ht.fn <- function (layer.prefix,layer.name,nonmissing) {
    # use this to make a quick plotting function
    layer <- raster(paste(layer.prefix,layer.name,sep=''))
    values(layer)[-nonmissing] <- NA # NOTE '-' NOT '!'
    load("../tort.coords.rasterGCS.Robj")
    ph <- function (x,...) { 
        values(layer)[nonmissing] <- x
        opar <- par()  # plotting layers messes up margins
        plot(layer,...)
        points(tort.coords.rasterGCS,pch=20,cex=.25)
        # par(opar[setdiff(names(opar), c("cin", "cra", "csi", "cxy", "din", "page") )])
        par(mar=opar$mar)
    }
    environment(ph) <- new.env()
    assign("tort.coords.rasterGCS",tort.coords.rasterGCS,environment(ph))
    return(ph)
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
