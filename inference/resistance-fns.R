#!/usr/bin/Rscript
require(Matrix)

# # find what directory this file is in
frame_files <- lapply(sys.frames(), function(x) x$ofile)
frame_files <- Filter(Negate(is.null), frame_files)
.PATH <- dirname(frame_files[[length(frame_files)]])
source(file.path(.PATH,"objective-functions.R"))
source(file.path(.PATH,"input-output-fns.R"))

# number of cores for parallel
getcores <- function (subdir) {
    # check in environment variables first
    if ( "parallel" %in% .packages()) {
        numcores.env <- as.numeric( Sys.getenv("MC_CORES",unset=NA) )
        if ( is.numeric(numcores.env) && !is.na(numcores.env) ) { 
            numcores <- numcores.env
        } else {
            numcores <- detectCores()
        }
    } else {
        numcores <- 1
    }
    if ( length(numcores)==0 || ( !missing(subdir) && ( as.numeric(gsub("x","",subdir)) < 50 ) ) ) {
        numcores <- 1
    }
    return(numcores)
}

###
# low-level functions for constructing e.g. adjacency matrices

# move @p to @j in dgCMatrix
p.to.j <- function (p) { rep( seq.int( length(p)-1 ), diff(p) ) }

make.G <- function (layer,nonmissing) {
    # set up the generator (G) matrix:
    #  layer is a raster layer 
    # and nonmissing is the indices of cells we want to be able to move to/from
    ij <- adjacent(layer,cells=nonmissing,target=nonmissing,directions=4,pairs=TRUE,sorted=TRUE) # to and from cells both loc
    ij <- ij[,2:1]
    stopifnot( all(ij[,1] != ij[,2]) ) ## NO DIAGONAL
    return( sparseMatrix( i=match(ij[,1],nonmissing), j=match(ij[,2],nonmissing), x=1.0 ) )
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

hitting.colinearity <- function (params, locs, obs.locs, G, update.G, layers, transfn, valfn, ndelta, ngamma, 
    hs=hitting.sensitivity(params, locs, G, update.G, layers, transfn, valfn, ndelta, ngamma, numcores=numcores),
    numcores=getcores()) {
    # use hitting.sensitivity to construct the Gram matrix
    # with rows and columns indexed by parameters
    # giving the amount of colinearity between the effects of each parameter on the hitting times at the observed locations
    obs.hs <- hs$gradient[obs.locs,,]
    dim(obs.hs) <- c( prod(dim(obs.hs)[1:2]), dim(obs.hs)[3] )
    return( cov(obs.hs) )
}

hitting.sensitivity <- function (params, neighborhoods, G, update.G, layers, transfn, valfn, ndelta, ngamma, do.hessian=FALSE, numcores=getcores()) {
    # return a list whose [[k]]th entry is the matrix of derivatives 
    #   of the hitting times with respect to the k-th parameter
    # and optionally the second derivatives also
    # FOR logistic transform
    stopifnot( transfn(2) == 1/(1+exp(-2)) )
    if ( numcores>1 && "parallel" %in% .packages()) {
        this.apply <- function (...) { do.call( cbind, mclapply( ..., mc.cores=numcores ) ) }
    } else {
        this.apply <- function (...) { sapply( ... ) }
    }
    nparams <- length(params)
    gamma <- params[1+(1:ngamma)]
    delta <- params[1 + ngamma + (1:ndelta)]
    G@x <- update.G(params)
    G.d <- G - Diagonal(nrow(G),rowSums(G))
    hts <- hitting.analytic(neighborhoods,G.d,numcores)
    zeros <- unlist(neighborhoods) + rep((seq_along(neighborhoods)-1)*nrow(hts),sapply(neighborhoods,length))
    GH <- G.d %*% hts
    GH[zeros] <- 0
    gradient <- numeric( nrow(G)*length(neighborhoods)*length(params) )
    dim(gradient) <- c( nrow(G), length(neighborhoods), length(params) )
    dimnames(gradient) <- list( NULL, names(neighborhoods), names(params) )
    # beta.grad <- list( .hsolve( neighborhoods, G.d, (-1)*GH ) )
    # names(beta.grad) <- names(params)[1]
    ## gradient[,,1] <- .hsolve( neighborhoods, G.d, (-1)*GH )  # duh:
    gradient[,,1] <- (-1)*hts
    # gamma.grad <- lapply( seq(1,length.out=ngamma), function (kg) {
    gfac <- (1-transfn(valfn(gamma)))  # from first deriv of logistic
    gfac2 <- gfac * (1-2*transfn(valfn(gamma)))  # from second deriv of logistic
    gamma.G <- lapply( seq(1,length.out=ngamma), function (kg) {
            (layers[,kg] * gfac) * GH
        } )
    for (kg in seq(1,length.out=ngamma)) {
        gradient[,,1+kg] <- .hsolve( neighborhoods, G.d, (-1)*gamma.G[[kg]], numcores=numcores )
    }
    # names(gamma.grad) <- names(params)[seq(2,length.out=ngamma)]
    delfac <- (1-transfn(valfn(delta)[G@i+1L]+valfn(delta)[Gjj]))
    delfac2 <- delfac * (1-2*transfn(valfn(delta)[G@i+1L]+valfn(delta)[Gjj]))
    delta.G <- lapply( seq(1,length.out=ndelta), function (kd) {
            GL <- G
            GL@x <- G@x * ( layers[Gjj,kd] + layers[G@i+1L,kd] ) * delfac
            return( GL - Diagonal(nrow(GL),rowSums(GL)) )
        } )
    # delta.grad <- lapply( seq(1,length.out=ndelta), function (kd) {
    delta.GL <- lapply( seq(1,length.out=ndelta), function (kd) {
            delta.G[[kd]] %*% hts
        } )
    for (kd in seq(1,length.out=ndelta)) {
        gradient[,,1+ngamma+kd] <- .hsolve( neighborhoods, G.d, (-1)*delta.GL[[kd]], numcores=numcores )
    }
    # names(delta.grad) <- names(params)[seq(2+ngamma,length.out=ndelta)]
    if (do.hessian) {
        hessian <- numeric( nrow(G) * length(neighborhoods) * nparams^2 )
        dim(hessian) <- c(nrow(G),length(neighborhoods),nparams,nparams)
        dimnames(hessian) <- list( NULL, names(neighborhoods), names(params), names(params) )
        # second deriv wrt beta : note that (d/d beta) hts = -hts
        ## hessian[,,1,1] <- (-1) * .hsolve( neighborhoods, G.d, 2 * (G.d %*% gradient[,,1]) + GH ) # duh:
        hessian[,,1,1] <- hts
        # beta.gamma.curv <- lapply( seq(1,length.out=ngamma), function (kg) {
        for (kg in seq(1,length.out=ngamma)) {
            hessian[,,1,1+kg] <- hessian[,,1+kg,1] <- (-1) * gradient[,,1+kg]
        }
        # beta.delta.curv <- lapply( seq(1,length.out=ndelta), function (kd) {
        for (kd in seq(1,length.out=ndelta)) {
            hessian[,,1,1+ngamma+kd] <- hessian[,,1+ngamma+kd,1] <- (-1) * gradient[,,1+ngamma+kd]
        }
        # gamma.curv <- lapply( seq(1,length.out=ngamma), function (kg1) {
        #     lapply( seq(1,length.out=ngamma), function (kg2) {
        # Note numerous simplifications below since (d/d gamma) G = (diagonal matrix) * G
        for (kg1 in seq(1,length.out=ngamma)) {
            hessian[,,1+kg1,1+kg1] <- .hsolve( neighborhoods, G.d,
                    ( 2 * layers[,kg1] * gfac) * gamma.G[[kg1]]
                    + (layers[,kg1]^2 * gfac2) 
                , numcores=numcores)
            for (kg2 in seq(1,length.out=kg1-1)) {
                hessian[,,1+kg1,1+kg2] <- hessian[,,1+kg2,1+kg1] <- .hsolve( neighborhoods, G.d,
                        (layers[,kg1] * gfac) * gamma.G[[kg2]]
                        + (layers[,kg2] * gfac) * gamma.G[[kg1]]
                        + (layers[,kg1] * layers[,kg2] * gfac2) 
                    , numcores=numcores)
            }
        }
        # delta.curv <- lapply( seq(1,length.out=ndelta), function (kd1) {
        #     lapply( seq(1,length.out=ndelta), function (kd2) {
        for (kd1 in seq(1,length.out=ndelta)) {
            GL2 <- G
            GL2@x <- G@x * ( layers[Gjj,kd1] + layers[G@i+1L,kd1] )^2 * delfac2
            dG2L <- rowSums(GL2)
            hessian[,,1+ngamma+kd1,1+ngamma+kd1] <- (-1) * .hsolve( neighborhoods, G.d,
                    2 * ( delta.G[[kd1]] %*% gradient[,,1+ngamma+kd1] )
                    + GL2 %*% hts - dG2L * hts
                , numcores=numcores)
            for (kd2 in seq(1,length.out=kd1-1)) {
                GL2 <- G
                GL2@x <- G@x * ( layers[Gjj,kd1] + layers[G@i+1L,kd1] ) * ( layers[Gjj,kd2] + layers[G@i+1L,kd2] ) * delfac2
                dG2L <- rowSums(GL2)
                hessian[,,1+ngamma+kd1,1+ngamma+kd2] <- hessian[,,1+ngamma+kd2,1+ngamma+kd1] <- (-1) * .hsolve( neighborhoods, G.d,
                        ( delta.G[[kd1]] %*% gradient[,,1+ngamma+kd2] )
                        + ( delta.G[[kd2]] %*% gradient[,,1+ngamma+kd1] )
                        + GL2 %*% hts - dG2L * hts
                    , numcores=numcores)
            }
        }
        # gamma.delta.curv <- lapply( seq(1,length.out=ngamma), function (kg) {
        #         lapply( seq(1,length.out=ndelta), function (kd) {
        for (kg in seq(1,length.out=ngamma)) {
            for (kd in seq(1,length.out=ndelta)) {
                hessian[,,1+kg,1+ngamma+kd] <- hessian[,,1+ngamma+kd,1+kg] <- (-1) * .hsolve( neighborhoods, G.d,
                            (-1) * (layers[,kg] * gfac) * delta.GL[[kd]]
                            + ( delta.G[[kd]] %*% gradient[,,1+kg] )
                            + (layers[,kg] * gfac) * ( delta.G[[kd]] %*% hts )
                        , numcores=numcores)
            }
        }
    } else { hessian <- NULL }
    return( list( gradient=gradient, hessian=hessian ) )
}

.hsolve <- function ( neighborhoods, G.d, x, numcores=1 ) {
    this.apply <- if ( numcores>1 && "parallel" %in% .packages()) {
                function (...) { do.call( cbind, mclapply( ..., mc.cores=numcores ) ) }
            } else { function (...) { sapply( ... ) } }
    this.apply( seq_along(neighborhoods), function (loc.ind) {
            klocs <- neighborhoods[[loc.ind]][!is.na(neighborhoods[[loc.ind]])]
            if (length(klocs)>0) {
                z <- numeric(nrow(G.d))
                z[-klocs] <- as.vector( solve( G.d[-klocs,-klocs], x[-klocs,loc.ind] ) )
                return( z )
            } else {
                return(NA)
            }
        } )
}

interp.hitting <- function ( locs, G, obs.ht, obs.locs, alpha=1, blocked=numeric(0), numcores=getcores() ) {
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
    if (is.list(blocked)) { blocked <- unlist(blocked) }
    if (length(blocked)>0) {
        Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )
        G@x[Gjj%in%blocked] <- 0
    }
    G <- G - Diagonal(nrow(G),rowSums(G))
    # Pmat projects full hitting times onto the obs.locs
    Pmat <- sparseMatrix( i=seq_along(obs.locs), j=obs.locs, x=1, dims=c(length(obs.locs),nrow(G)) )
    PtP <- alpha * crossprod(Pmat)
    hts <- this.apply( seq_along(locs), function (k) {
            klocs <- unique(c( unlist(locs[k])[!is.na(unlist(locs[k]))], blocked ))
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

J.from.G <- function (params,layers,transfn,G) {
    pivec <- stationary.dist( params, layers, transfn )
    # try with the symmetrized matrix
    Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )
    J <- G
    J@x <- G@x * ( sqrt(pivec)[1L+G@i] / sqrt(pivec)[Gjj] )
    J <- forceSymmetric(J)
    return( list(pivec=pivec, J=J) )
}

interp.hitting.sym <- function ( locs, J, pivec, obs.ht, obs.locs, alpha=1, blocked=numeric(0), numcores=getcores() ) {
    # as interp.hitting
    #   but G = pivec^(-1/2) J pivec^(1/2)
    if ( numcores>1 && "parallel" %in% .packages()) {
        this.apply <- function (...) { do.call( cbind, mclapply( ..., mc.cores=numcores ) ) }
    } else {
        this.apply <- function (...) { sapply( ... ) }
    }
    if (is.list(blocked)) { blocked <- unlist(blocked) }
    if (length(blocked)>0) {
        Jjj <- rep( seq.int(length(J@p)-1), diff(J@p) )
        J@x[Jjj%in%blocked] <- 0
    }
    sqrt.pi <- sqrt(pivec)
    nu <- 1/sqrt.pi
    Jd <- J - Diagonal( nrow(J), as.vector(sqrt.pi*(J%*%nu)) )
    obs.nu <- obs.ht*nu[obs.locs]
    # Pmat projects full hitting times onto the obs.locs
    Pmat <- sparseMatrix( i=seq_along(obs.locs), j=obs.locs, x=1, dims=c(length(obs.locs),nrow(Jd)) )
    PtP <- alpha * crossprod(Pmat)
    hts <- this.apply( seq_along(locs), function (k) {
            klocs <- unique(c( unlist(locs[k])[!is.na(unlist(locs[k]))], blocked ))
            if (length(klocs)>0) {
                bvec <- as.vector( alpha * crossprod(Pmat[,-klocs],obs.nu[,k]) + crossprod( Jd[-klocs,-klocs], (-1)*nu[-klocs] ) )
                z <- numeric(nrow(G))
                z[-klocs] <- sqrt.pi[-klocs] * as.vector( solve( PtP[-klocs,-klocs] + crossprod(Jd[-klocs,-klocs]), bvec ) )
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
    this.lapply <- if ( numcores>1 && "parallel" %in% .packages()) { function (...) { mclapply( ..., mc.cores=numcores ) } } else { lapply }
    this.lapply( seq_along(neighborhoods), function (k) {
            nh <- neighborhoods[[k]]
            bd <- boundaries[[k]]
            as.matrix( solve( G[nh,nh]-Diagonal(n=length(nh),x=dG[nh]), -G[nh,bd,drop=FALSE] ) )
        } )
}

get.hitting.times <- function (G,dG,neighborhoods,boundaries,numcores=getcores()) {
    # returns a list of vectors with the [[k]]th has [i]th entry
    # the hitting times from the i-th element of neighborhoods[[k]] to boundaries[[k]]
    #  (like hitting.analytic but different syntax)
    this.lapply <- if ( numcores>1 && "parallel" %in% .packages()) { function (...) { mclapply( ..., mc.cores=numcores ) } } else { lapply }
    this.lapply( seq_along(neighborhoods), function (k) {
            nh <- neighborhoods[[k]]
            bd <- boundaries[[k]]
            as.vector( solve( G[nh,nh]-Diagonal(n=length(nh),x=dG[nh]), rep.int(-1.0,length(nh)) ) )
        } )
}


####
# model things

stationary.dist <- function (params,layers,transfn) {
    # proportional to the stationary distribution
    gamma <- params[2:(1+ncol(layers))]
    return( 1 / transfn( rowSums( layers * gamma[col(layers)] ) ) )
}

########
# Raster whatnot

get.neighborhoods <- function ( ndist, locations, nonmissing, layer, numcores=getcores(), na.rm=TRUE ) {
    # locations should be either a SpatialPoints object or a 2-column matrix of coordinates
    this.lapply <- if ( numcores>1 && "parallel" %in% .packages()) { function (...) { mclapply( ..., mc.cores=numcores ) } } else { lapply }
    if ( class(locations)=="SpatialPoints" ) { locations <- coordinates(locations) }  # this is the first thing distanceFromPoints does anyhow
    if (is.null(dim(locations))) { locations <- matrix(locations,ncol=2) }
    neighborhoods <- this.lapply( 1:NROW(locations) , function (k) {
            d_tort <- distanceFromPoints( layer, locations[k,] ) 
            match( Which( d_tort <= max(ndist,minValue(d_tort)), cells=TRUE, na.rm=TRUE ), nonmissing )
        } )
    if (na.rm) { neighborhoods <- lapply(neighborhoods,function (x) { x[!is.na(x)] }) }
    return(neighborhoods)
}


get.boundaries <- function ( neighborhoods, nonmissing, layer, numcores=getcores(), na.rm=TRUE ) {
    this.lapply <- if ( numcores>1 && "parallel" %in% .packages()) { function (...) { mclapply( ..., mc.cores=numcores ) } } else { lapply }
    boundaries <- this.lapply( neighborhoods, function (nh) {
            values(layer) <- TRUE
            values(layer)[nonmissing][nh] <- NA
            bdry <- boundaries(layer,directions=4)
            match( which( (!is.na(values(bdry))) & (values(bdry)==1) ), nonmissing )
        } )
    if (na.rm) { boundaries <- lapply(boundaries,function (x) { x[!is.na(x)] }) }
    return(boundaries)
}

remove.clumps <- function (layer) {
    # NA out everything but the biggest connected cluster
    cl <- clump(layer,directions=4)
    cl.table <- table(values(cl))
    big.clump <- as.numeric(names(cl.table))[ which.max( table( values(cl) ) ) ]
    layer[cl!=big.clump] <- NA
    return(layer)
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
