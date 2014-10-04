#!/usr/bin/Rscript

###
# Get hitting times with a landscape layer
#   e.g.
#     Rscript initial-hitting-times.R ../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_ annual_precip
# (note the space!)                                                                           ---->^

source("resistance-fns.R")
require(raster)

require(parallel)
numcores<-as.numeric(scan(pipe("cat /proc/cpuinfo | grep processor | tail -n 1 | awk '{print $3}'")))+1

if (!interactive()) {
    layer.prefix <- commandArgs(TRUE)[1]
    subdir <- commandArgs(TRUE)[2]
    layer.file <- commandArgs(TRUE)[3]
    param.file <- commandArgs(TRUE)[4] 
    method <- commandArgs(TRUE)[5] 
    prev.ht <- if (length(commandArgs(TRUE))>5) { commandArgs(TRUE)[6] } else { NULL } 
    maxtime <- if (length(commandArgs(TRUE))>6) { commandArgs(TRUE)[7] } else { 6*60*60 } 
    outfile <- if (length(commandArgs(TRUE))>7) { commandArgs(TRUE)[8] } else { NULL }
} else {
    layer.prefix <- "../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_"
    subdir <- "100x"
    layer.file <- "../inference/six-raster-list"
    param.file <- "simple-init-params-six-raster-list.tsv"
    method <- "CG"
    prev.ht <- "100x/500x-aggregated-hitting-times.tsv"
    maxtime <- 2*60
    outfile <- NULL

    # layer.prefix <- "../geolayers/TIFF/500x/500x_"
    # subdir <- "500x"
    # layer.file <- "../inference/six-raster-list"
    # param.file <- "simple-init-params-six-raster-list.tsv"
    # method <- "analytic"
    # prev.ht <- NULL
    # maxtime <- 1*60
    # outfile <- NULL
}

if (is.null(outfile)) { paste( subdir, "/", basename(layer.file), "-hitting-times.tsv", sep='') }

layer.names <- scan(layer.file,what="char") 

load( paste(subdir,"/",basename(layer.prefix),"_",basename(layer.file),"_","G.RData",sep='') ) # provides "G"        "update.G" "ndelta"   "ngamma"   "transfn"  "valfn"    "layers"
Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )

load(paste(subdir,"/",basename(layer.prefix),"tortlocs.RData",sep='')) # provides 'locs'
na.indiv <- which( is.na( locs ) )
locs <- locs[-na.indiv]

##
# initial parameters?
#   in time t, 1D RW with rate r does t*r jumps, displacement has variance t*r
#   so time to move N grid sites away is sqrt(N)/r
#   so if hitting times are of order T, want r of order sqrt(N)/T

init.param.table <- read.table( param.file, header=TRUE )
init.params <- unlist( init.param.table[ match( subdir, init.param.table[,1] ), -1 ] )

G@x <- update.G(init.params)

if (method=="analytic") {

    hts <- hitting.analytic( locs, G-diag(rowSums(G)), numcores=numcores )

} else if (method=="CG") {

    init.hts <- as.matrix( read.table(prev.ht,header=TRUE) )

    stopifnot( nrow(init.hts) == nrow(G) )
    # solve for hitting times

    dG <- rowSums(G)
    cG <- colSums(G)
    # objective function
    H <- function (ht,loc) {
        # ( (G-D)ht + 1 )^T S ( (G-D)ht + 1)
        # where S = I except S[loc,loc]=0
        ht[loc] <- 0
        z <- G%*%ht - dG*ht + 1
        z[loc] <- 0
        return( ( sum( z^2 ) )/length(z) )
    }
    dH <- function (ht,loc) {
        # 2 (G-D)^T S ( (G-D) ht + 1 )
        ht[loc] <- 0
        z <- G%*%ht - dG*ht + 1
        z[loc] <- 0
        z <- (crossprod(G,z) - dG*z)
        z[loc] <- 0
        return( 2 * as.vector(z) / length(z) )
    }

    # parscale <- rep( nrow(G) / exp( mean( log(dG), trim=.1, na.rm=TRUE ) ), nrow(G) )
    parscale <- rep( mean(init.hts), nrow(init.hts) )

    H.time <- system.time( sapply( 1:ncol(init.hts), function (k) H(init.hts[,k],loc=locs[k]) ) ) / ncol(init.hts)
    dH.time <- system.time( sapply( 1:ncol(init.hts), function (k) dH(init.hts[,k],loc=locs[k]) ) ) / ncol(init.hts)

    maxit <- floor( maxtime / (H.time[1] + dH.time[1]) * numcores / ncol(init.hts) ) / 10  # turns out optim adds in a fair bit of overhead, hence the '/10'

    optim.ht.list <- mclapply( seq_along(locs), function (loc.ind) {
                optim( par=init.hts[,loc.ind], fn=H, gr=dH, loc=locs[loc.ind], 
                    method="L-BFGS-B", control=list( parscale=parscale, maxit=maxit ), lower=0, upper=Inf ) 
            }, mc.cores=numcores )

    convergences <- sapply(optim.ht.list,"[[","convergence")
    unconverged <- which(convergences != 0)

    save( optim.ht.list, file=gsub(".tsv", "-optim.RData", outfile) )

    hts <- sapply( optim.ht.list, "[[", "par" )

    if (FALSE) {
        # check gradient
        loc.ind <- 10
        loc <- locs[loc.ind]
        ht <- init.hts[,loc.ind]

        eps <- 1e-5 * runif(length(ht))
        c( H(ht,loc=loc), H(ht+eps,loc=loc)-H(ht,loc=loc), sum(eps*dH(ht,loc=loc)) )

        # look at convergence
        load( paste(subdir, "/", basename(layer.prefix),"_", basename(layer.file),"_nonmissing.RData",sep='') ) # provides nonmissing
        ph <- plot.ht.fn(layer.prefix,"annual_precip",nonmissing)

        loc.ind <- 10
        oht.list <- vector( mode='list', length=6 )
        oht.list[[1]] <- list( par=init.hts[,loc.ind] )
        for (k in 2:10) { oht.list[[k]] <- optim( par=oht.list[[k-1]]$par, fn=H, gr=dH, loc=locs[loc.ind], method="L-BFGS-B", control=list( parscale=parscale, maxit=100 ), lower=0, upper=Inf )  }

        layout( matrix(1:6,nrow=2) )
        for (k in 2:length(oht.list)) { 
            ph( oht.list[[k]]$par )
            ph( oht.list[[k]]$par - oht.list[[k-1]]$par )
            ph( dH( oht.list[[k-1]]$par, loc=locs[loc.ind] ) )
            if (is.null(locator(1))) break
        }

        layout( matrix(1:6,nrow=2) )
        ph( dH( oht.list[[length(oht.list)]]$par, loc=locs[loc.ind] ) )
        for (k in floor(seq(1,length(oht.list)-1,length.out=5))) {
            ph( oht.list[[length(oht.list)]]$par - oht.list[[k]]$par, main=k )
        }

        # look at inferred layer
        ph( (transfn(valfn(init.params[1 + (1:ngamma)]))) )
        ph( log(transfn(valfn(init.params[1 + (1:ngamma)]))) )
        ph( (transfn(valfn(init.params[1 + ngamma + (1:ndelta)]))) )
        ph( log(transfn(valfn(init.params[1 + ngamma + (1:ndelta)]))) )

        # compare answers to initial values
        layout( matrix(1:6,nrow=2) )
        par(mar=c(5,4,4,5)+.1)
        for (k in 1:ncol(hts)) {
            diff.hts <- hts[,k] - init.hts[,k]
            diff.hts[ (diff.hts<quantile(diff.hts,.05,na.rm=TRUE)) | (diff.hts>quantile(diff.hts,.95,na.rm=TRUE)) ] <- NA
            ph( init.hts[,k] )
            ph( hts[,k] )
            ph( diff.hts )
            if (is.null(locator(1))) break
        }

    }
}

colnames( hts ) <- locs
write.table( hts, file=outfile, row.names=FALSE )

if (FALSE) {

    load( paste(subdir, "/", basename(layer.prefix),"_", basename(layer.file),"_nonmissing.RData",sep='') ) # provides nonmissing
    ph <- plot.ht.fn(layer.prefix,"annual_precip",nonmissing)

    layout( matrix(1:6,nrow=2) )
    for (k in 1:ncol(hts)) { ph( hts[,k] ); if ((k%%6==0) && is.null(locator(1))) break }

}

if (FALSE) {

    ####
    ## FIRST
    # infer parameters that better match these hitting times

    # Massage the numerics.
    zeros <- which( row(init.hts) == locs[col(init.hts)] )
    scaling <- sqrt(nrow(G) * length(locs))
    init.hts <- init.hts/scaling
    init.hts[zeros] <- 0
    sc.one <- 1/scaling

    dG <- rowSums(G)
    GH <- G %*% init.hts - dG*init.hts
    GH[zeros] <- 0

    init.beta <- init.params[1]

    L <- function (params) {
        if (any(params != get("params",parent.env(environment()) ) ) ) { 
            assign("params", params, parent.env(environment()) )
            evalq( G@x <- update.G(c(init.beta,params)), parent.env(environment()) )
            evalq( dG <- rowSums(G), parent.env(environment()) )
            GH <- G %*% init.hts - dG*init.hts
            GH[zeros] <- 0
            assign("GH", GH, parent.env(environment()) )
        }
        ans <- ( sum( (GH+sc.one)^2 ) - length(zeros)*sc.one^2 )
        if (!is.finite(ans)) { browser() }
        return(ans)
    }
    dL <- function (params) {
        if (any(params != get("params", parent.env(environment()) ) ) ) { 
            assign("params", params, parent.env(environment()) )
            evalq( G@x <- update.G(c(init.beta,params)), parent.env(environment()) )
            evalq( dG <- rowSums(G), parent.env(environment()) )
            GH <- G %*% init.hts - dG*init.hts
            GH[zeros] <- 0
            assign("GH", GH, parent.env(environment()) )
        }
        ggrads <- sapply( 1:ncol(layers), function (kk) {
                2 * sum( layers[,kk] * GH * (GH+sc.one) )
            } )
        dgrads <- ggrads + sapply( 1:ncol(layers), function (kk) {
                GL <- G
                GL@x <- G@x * layers[Gjj,kk]
                dGL <- rowSums(GL)
                GLH <- GL %*% init.hts - dGL*init.hts
                GLH[zeros] <- 0
                return( 2 * sum( GLH * (GH+sc.one)  ) )
            } )
        ans <- ( c(ggrads, dgrads) )
        if (any(!is.finite(ans))) { browser() }
        return(ans)
    }
    environment(L) <- environment(dL) <- fun.env <- list2env( list(
                    G=G,
                    dG=dG,
                    params=init.params[-1],
                    GH=GH), 
            parent=environment() )

    L(init.params[-1])
    dL(init.params[-1])

    L(init.params[-1]+.01)
    dL(init.params[-1]+.01)

    parscale <- c( rep(0.1,length(init.params)-1) )
    results <- optim( par=init.params[-1], fn=L, gr=dL, control=list(parscale=parscale), method="BFGS" )
    if (results$convergence != 0) {
        results <- optim( par=results$par, fn=L, gr=dL, control=list(parscale=parscale/10), method="BFGS" )
    }


    ####
    ## THEN

}
