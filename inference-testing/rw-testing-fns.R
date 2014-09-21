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

