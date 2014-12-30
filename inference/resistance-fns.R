#!/usr/bin/Rscript

# # find what directory this file is in
frame_files <- lapply(sys.frames(), function(x) x$ofile)
frame_files <- Filter(Negate(is.null), frame_files)
.PATH <- dirname(frame_files[[length(frame_files)]])
source(file.path(.PATH,"objective-functions.R"))
source(file.path(.PATH,"input-output-fns.R"))

# number of cores for parallel
getcores <- function (subdir) {
    if ( "parallel" %in% .packages()) {
        cpupipe <- pipe("cat /proc/cpuinfo | grep processor | tail -n 1 | awk '{print $3}'")
        numcores <- 1+as.numeric(scan(cpupipe))
        close(cpupipe)
    } else {
        numcores <- 1
    }
    if ( !missing(subdir) && ( as.numeric(gsub("x","",subdir)) < 50 ) ) {
        numcores <- 1
    }
    return(numcores)
}

###
# low-level functions for constructing e.g. adjacency matrices

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


###
#  Functions for solving hitting times and the like

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


hitting.analytic <- function (locs, G, numcores=getcores(), blocked=numeric(0)) {
    # compute analytical expected hitting times
    #   here `locs` is a vector of (single) locations
    #     or a list of vectors
    #   G is a generator matrix WITHOUT diagonal
    # Optionally, block all movement through celled indexed by 'blocked'
    if ( numcores>1 && "parallel" %in% .packages()) {
        this.apply <- function (...) { do.call( cbind, mclapply( ..., mc.cores=numcores ) ) }
    } else {
        this.apply <- function (...) { sapply( ... ) }
    }
    # zero out entries into or out of cells indexed in 'blocked'
    if (is.list(blocked)) { blocked <- unlist(blocked) }
    if (length(blocked)>0) {
        Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )
        G@x[Gjj%in%blocked] <- 0
    }
    G <- G - Diagonal(nrow(G),rowSums(G))
    hts <- this.apply( locs, function (k) { 
                klocs <- unique(c( k[!is.na(k)], blocked ))
                if (length(klocs)>0) {
                    z <- numeric(nrow(G))
                    z[-klocs] <- as.vector( solve( G[-klocs,-klocs], rep.int(-1.0,nrow(G)-length(klocs)) ) )
                    return( z )
                } else {
                    return(NA)
                }
            } )
    return(hts)
}

hitting.colinearity <- function (params, locs, obs.locs, G, update.G, layers, transfn, valfn, ndelta, ngamma, numcores=getcores()) {
    # use hitting.sensitivity to construct the Gram matrix
    # with rows and columns indexed by parameters
    # giving the amount of colinearity between the effects of each parameter on the hitting times at the observed locations
    hs <- hitting.sensitivity(params, locs, G, update.G, layers, transfn, valfn, ndelta, ngamma, numcores=numcores)
    obs.hs <- sapply( hs, function (x) { x[obs.locs,] } )
    return( cov(obs.hs) )
}

hitting.sensitivity <- function (params, locs, G, update.G, layers, transfn, valfn, ndelta, ngamma, numcores=getcores()) {
    # return a list whose [[k]]th entry is the matrix of derivatives 
    #   of the hitting times with respect to the k-th parameter
    # FOR logistic transform
    stopifnot( transfn(2) == 1/(1+exp(-2)) )
    if ( numcores>1 && "parallel" %in% .packages()) {
        this.apply <- function (...) { do.call( cbind, mclapply( ..., mc.cores=numcores ) ) }
    } else {
        this.apply <- function (...) { sapply( ... ) }
    }
    gamma <- params[1+(1:ngamma)]
    delta <- params[1 + ngamma + (1:ndelta)]
    G@x <- update.G(params)
    G.d <- G - Diagonal(nrow(G),rowSums(G))
    hts <- hitting.analytic(locs,G.d,numcores)
    zeros <- unlist(locs) + rep((seq_along(locs)-1)*nrow(hts),sapply(locs,length))
    GH <- G.d %*% hts
    GH[zeros] <- 0
    bgrad <- list( this.apply( seq_along(locs), function (loc.ind) {
                klocs <- locs[[loc.ind]][!is.na(locs[[loc.ind]])]
                if (length(klocs)>0) {
                    z <- numeric(nrow(G.d))
                    z[-klocs] <- (-1) * as.vector( solve( G.d[-klocs,-klocs], GH[-klocs,loc.ind] ) )
                    return( z )
                } else {
                    return(NA)
                }
            } ) )
    names(bgrad) <- names(params)[1]
    ggrad <- lapply( seq(1,length.out=ngamma), function (kg) {
            LGH <- (layers[,kg] * (1-transfn(valfn(gamma)))) * GH
            this.apply( seq_along(locs), function (loc.ind) {
                    klocs <- locs[[loc.ind]][!is.na(locs[[loc.ind]])]
                    if (length(klocs)>0) {
                        z <- numeric(nrow(G.d))
                        z[-klocs] <- (-1) * as.vector( solve( G.d[-klocs,-klocs], LGH[-klocs,loc.ind] ) )
                        return( z )
                    } else {
                        return(NA)
                    }
                } )
        } )
    names(ggrad) <- names(params)[seq(2,length.out=ngamma)]
    dgrad <- lapply( seq(1,length.out=ndelta), function (kd) {
            GL <- G
            GL@x <- G@x * ( layers[Gjj,kd] + layers[G@i+1L,kd] ) * (1-transfn(valfn(delta)[G@i+1L]+valfn(delta)[Gjj]))
            dGL <- rowSums(GL)
            GLH <- GL %*% hts - dGL*hts
            GLH[zeros] <- 0
            this.apply( seq_along(locs), function (loc.ind) {
                    klocs <- locs[[loc.ind]][!is.na(locs[[loc.ind]])]
                    if (length(klocs)>0) {
                        z <- numeric(nrow(G))
                        z[-klocs] <- (-1) * as.vector( solve( G.d[-klocs,-klocs], GLH[-klocs,loc.ind] ) )
                        return( z )
                    } else {
                        return(NA)
                    }
                } )
        } )
    names(dgrad) <- names(params)[seq(2+ngamma,length.out=ndelta)]
    return( c( bgrad, ggrad, dgrad ) )
}

interp.hitting <- function ( locs, G, obs.ht, obs.locs, alpha=1, numcores=getcores() ) {
    # interpolate hitting times by minimizing squared error:
    #       G is a generator matrix WITHOUT diagonal
    #       locs is a vector of indices, or a list of vectors, of the rows of G for which we have data
    #       obs.ht is the (locs x locs) matrix of mean hitting times
    #       obs.locs is the indices for which obs.ht correspond
    #       alpha is a fudge factor (so far unnecessary?)
    #   note that without the second 'crossprod' term in the defn of 'bvec'
    #   this finds the harmonic function interpolating the observed hitting times
    #     (or something close, I think)
    if ( numcores>1 && "parallel" %in% .packages()) {
        this.apply <- function (...) { do.call( cbind, mclapply( ..., mc.cores=numcores ) ) }
    } else {
        this.apply <- function (...) { sapply( ... ) }
    }
    G <- G - Diagonal(nrow(G),rowSums(G))
    # Pmat projects full hitting times onto the obs.locs
    Pmat <- sparseMatrix( i=seq_along(obs.locs), j=obs.locs, x=1, dims=c(length(obs.locs),nrow(G)) )
    PtP <- alpha * crossprod(Pmat)
    hts <- this.apply( seq_along(locs), function (k) {
            klocs <- unlist(locs[k])[!is.na(k)]
            if (length(klocs)>0) {
                bvec <- as.vector( alpha * crossprod(Pmat[,-klocs],obs.ht[,k]) + crossprod( G[-klocs,-klocs], rep(-1.0,nrow(G)-length(klocs)) ) )
                z <- numeric(nrow(G))
                z[-klocs] <- as.vector( solve( PtP[-klocs,-klocs] + crossprod(G[-klocs,-klocs]), bvec ) )
                return( z )
            } else {
                return( NA )
            }
    } )
    return(hts)
}

interp.tradeoff <- function ( hts, locs, G, dG, obs.ht, obs.locs, numcores=getcores() ) {
    # Evaluate the two quantities that are being minimized by interp.hitting:
    #   | G %*% hts + 1 |  and  | hts[obs.locs,] - obs.ht |
    stopifnot( all.equal(diag(G),numeric(nrow(G))) )
    zeros <- unlist(locs) + rep((seq_along(locs)-1)*nrow(hts),sapply(locs,length))
    hts[zeros] <- 0
    GH <- G %*% hts - dG * hts
    GH[zeros] <- 0
    return( c( sqrt(sum( (GH+1)^2 ) - length(zeros) ), sqrt(sum( (hts[obs.locs,] - obs.ht)^2 ) ) ) )
}

get.hitting.probs <- function (G,dG,neighborhoods,boundaries,numcores=getcores()) {
    # returns a list of matrices of with the [[k]]th has [i,j]th entry
    # the hitting probabilty from the i-th element of neighborhoods[[k]] to the j-th element of boundaries[[k]]
    mclapply( seq_along(neighborhoods), function (k) {
            nh <- neighborhoods[[k]]
            bd <- boundaries[[k]]
            as.matrix( solve( G[nh,nh]-Diagonal(n=length(nh),x=dG[nh]), -G[nh,bd,drop=FALSE] ) )
        }, mc.cores=numcores )
}

get.hitting.times <- function (G,dG,neighborhoods,boundaries,numcores=getcores()) {
    # returns a list of vectors with the [[k]]th has [i]th entry
    # the hitting times from the i-th element of neighborhoods[[k]] to boundaries[[k]]
    #  (like hitting.analytic but different syntax)
    mclapply( seq_along(neighborhoods), function (k) {
            nh <- neighborhoods[[k]]
            bd <- boundaries[[k]]
            as.vector( solve( G[nh,nh]-Diagonal(n=length(nh),x=dG[nh]), rep.int(-1.0,length(nh)) ) )
        }, mc.cores=numcores )
}


########
# Raster whatnot

get.neighborhoods <- function ( ndist, locations, nonmissing, layer, numcores=getcores(), na.rm=TRUE ) {
    # locations should be either a SpatialPoints object or a 2-column matrix of coordinates
    if ( class(locations)=="SpatialPoints" ) { locations <- coordinates(locations) }  # this is the first thing distanceFromPoints does anyhow
    if (is.null(dim(locations))) { locations <- matrix(locations,ncol=2) }
    neighborhoods <- mclapply( 1:NROW(locations) , function (k) {
        d_tort <- distanceFromPoints( layer, locations[k,] ) 
        match( Which( d_tort <= max(ndist,minValue(d_tort)), cells=TRUE, na.rm=TRUE ), nonmissing )
    }, mc.cores=numcores )
    if (na.rm) { neighborhoods <- lapply(neighborhoods,function (x) { x[!is.na(x)] }) }
    return(neighborhoods)
}


get.boundaries <- function ( neighborhoods, nonmissing, layer, numcores=getcores(), na.rm=TRUE ) {
    boundaries <- mclapply( neighborhoods, function (nh) {
        values(layer) <- TRUE
        values(layer)[nonmissing][nh] <- NA
        bdry <- boundaries(layer,directions=4)
        match( which( (!is.na(values(bdry))) & (values(bdry)==1) ), nonmissing )
    }, mc.cores=numcores )
    if (na.rm) { boundaries <- lapply(boundaries,function (x) { x[!is.na(x)] }) }
    return(boundaries)
}

which.nonoverlapping <- function (neighborhoods) {
    # find a set of neighborhoods that are mutually nonoverlapping
    perm <- sample(length(neighborhoods))
    goodones <- rep(FALSE,length(perm))
    goodones[1] <- TRUE
    for (k in seq_along(perm)[-1]) {
        goodones[k] <- ( 0 == length( intersect( neighborhoods[[k]], unlist(neighborhoods[goodones]) ) ) )
    }
    return( which(goodones) )
}

LinesFromIndices <- function ( ind.pairs, locations, layer, proj4string=CRS(sp::proj4string(layer)) ) {
    # Returns a SpatialLines-class list of Lines objects, one for each pair of indices in ind.pairs,
    # with these indexing the poitns in sample.coords.
    #  ind.pairs : a two-column matrix of indices
    #  locations : a two-column matrix of coordinates
    #  layer : something to extract a proj4string from
    #  proj4string : a CRS object with a valid proj4string
    if ( class(locations)=="SpatialPoints" ) { locations <- coordinates(locations) }  # this is the first thing distanceFromPoints does anyhow
    SpatialLines( apply( ind.pairs, 1, function (ij) {
            Lines( Line( rbind( locations[ij[1],], locations[ij[2],] ) ), ID=paste(ij,collapse=":") )
        } ), proj4string=proj4string )
}

##
# move between resolutions

upsample <- function ( layer.vals, ag.fact, layer.1, nonmissing.1, layer.2, nonmissing.2, checkit=FALSE ) {
    # moves from layer.1 to layer.2, which must be related by a factor of ag.fact
    values(layer.1)[nonmissing.1] <- layer.vals
    layer.1.dis <- crop( disaggregate( layer.1, fact=ag.fact, method='bilinear' ), layer.2 )
    stopifnot( all( dim(layer.1.dis)==dim(layer.2) ) )
    # can skip this step, hopefully
    if (checkit) {
        layer.1.dis.res <- resample( layer.1.dis, layer.2 )
        stopifnot( all( abs( values(layer.1.dis)[nonmissing.2] - values(layer.1.dis.res)[nonmissing.2] ) < 1e-3 ) )
    }
    # get values out
    return( values(layer.1.dis)[nonmissing.2] )
}

upsample.hts <- function ( hts, ..., numcores=getcores() ) {
    new.hts <- do.call( cbind, mclapply( 1:ncol(hts), function (k) {
                upsample( hts[,k], ... )
        }, mc.cores=numcores ) )
    colnames(new.hts) <- colnames(hts)
    return(new.hts)
}

downsample <- function ( layer.vals, ag.fact, layer.1, nonmissing.1, layer.2, nonmissing.2, checkit=FALSE ) {
    # moves from layer.2 to layer.1, which must be related by a factor of ag.fact
    values(layer.2)[nonmissing.2] <- layer.vals
    layer.2.ag <- crop( aggregate( layer.2, fact=ag.fact, fun=mean, na.rm=TRUE ), layer.1 )
    stopifnot( all( dim(layer.2.ag)==dim(layer.1) ) )
    # get values out
    return( values(layer.2.ag)[nonmissing.1] )
}

downsample.hts <- function ( hts, ..., numcores=getcores() ) {
    new.hts <- do.call( cbind, mclapply( 1:ncol(hts), function (k) {
                downsample( hts[,k], ... )
        }, mc.cores=numcores ) )
    colnames(new.hts) <- colnames(hts)
    return(new.hts)
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
