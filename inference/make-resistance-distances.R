#!/usr/bin/Rscript

usage <- '
Get hitting times with a list of landscape layers:
        Rscript make-resistance-distances.R (layer prefix) (subdir) (layer file) (parameter file) (method) [initial guess] [max running time] [output file]
e.g.
        Rscript make-resistance-distances.R ../geolayers/multigrid/512x/crm_ 256x six-raster-list multigrid-six-raster-list.tsv numeric 256x/512x-six-raster-list-aggregated-hitting-times.tsv 120

    Here `method` is either "analytic" or "numeric".
'

if (!interactive()) {
    if (length(commandArgs(TRUE))<5) { stop(usage) }
    print(commandArgs(TRUE))
    layer.prefix <- commandArgs(TRUE)[1]
    subdir <- commandArgs(TRUE)[2]
    layer.file <- commandArgs(TRUE)[3]
    param.file <- commandArgs(TRUE)[4] 
    method <- commandArgs(TRUE)[5] 
    prev.ht <- if (length(commandArgs(TRUE))>5) { commandArgs(TRUE)[6] } else { NULL } 
    maxit <- if (length(commandArgs(TRUE))>6) { as.numeric(commandArgs(TRUE)[7]) } else { 100 } 
    outfile <- if (length(commandArgs(TRUE))>7) { commandArgs(TRUE)[8] } else { NULL }
} else {
    layer.prefix <- "../geolayers/multigrid/256x/crm_"
    subdir <- "256x"
    layer.file <- "dem-layer-list"
    param.file <- "params-dem-layer-list.tsv"
    method <- "analytic"
    prev.ht <- NULL
    maxit <- 100
    outfile <- NULL

    # layer.prefix <- "../geolayers/TIFF/500x/500x_"
    # subdir <- "500x"
    # layer.file <- "../inference/six-raster-list"
    # param.file <- "simple-init-params-six-raster-list.tsv"
    # method <- "analytic"
    # prev.ht <- NULL
    # maxit <- NULL
    # outfile <- NULL
    argvec <- scan(what='char')
    layer.prefix <- argvec[1]
    subdir <- argvec[2]
    layer.file <- argvec[3]
    param.file <- argvec[4]
    method <- argvec[5]
    prev.ht <- argvec[6]
    maxit <- argvec[7]
    outfile <- argvec[8]
}

# number of scaling & shifting steps
nscale <- 0

cat(paste(commandArgs()),"\n")

if ( length(commandArgs(TRUE))<5 ) { cat(usage); q() }
if (! method %in% c("analytic","numeric") ) { cat("  Method must be either 'analytic' or 'numeric'. \n"); stop(usage) }

source("resistance-fns.R")
require(raster)

require(parallel)
numcores <- getcores()

if (!exists("outfile")||is.null(outfile)||is.na(outfile)) { outfile <- paste( subdir, "/", basename(layer.file), "-hitting-times.tsv", sep='') }

layer.names <- scan(layer.file,what="char") 

load( paste(subdir,"/",basename(layer.prefix),"_",basename(layer.file),"_","G.RData",sep='') ) # provides "G"    "Gjj"    "update.G" "ndelta"   "ngamma"   "transfn"  "valfn"    "layers"

load( paste( subdir, "/", basename(layer.prefix), basename(layer.file), "_neighborhoods.RData", sep='' ) ) # provides 'neighborhoods'
load(paste(subdir,"/",basename(layer.prefix),"tortlocs.RData",sep='')) # provides 'locs'

# REMOVE MISSING INDIV
na.indiv <- which( is.na( locs ) )
locs <- locs[-na.indiv]
neighborhoods <- lapply(neighborhoods[-na.indiv],function (x) { x[!is.na(x)] })


##
# initial parameters?
#   in time t, 1D RW with rate r does t*r jumps, displacement has variance t*r
#   so time to move N grid sites away is sqrt(N)/r
#   so if hitting times are of order T, want r of order sqrt(N)/T

init.param.table <- read.table( param.file, header=TRUE )
init.params <- unlist( init.param.table[ match( subdir, init.param.table[,1] ), -1 ] )

G@x <- update.G(init.params)

if (method=="analytic") {

    hts <- hitting.analytic( neighborhoods, G-diag(rowSums(G)), numcores=numcores )

} else if (method=="numeric") {

    init.hts <- as.matrix( read.table(prev.ht,header=TRUE) )
    init.hts <- sapply( 1:ncol(init.hts), function (k) {
                x <- init.hts[,k]
                x[neighborhoods[[k]]] <- 0
                return(x) } )

    stopifnot( nrow(init.hts) == nrow(G) )

    ###
    # solve for hitting times

    dG <- rowSums(G)
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
    parscale <- rep(mean(init.hts),nrow(init.hts)) # (mean(init.hts)+init.hts)/2

    # First get the overall scaling right:
    #  (d/da) | a Ax - b |^2 = 2 x^T A^T ( a Ax - b )
    #    = 0  =>  a = x^T A b / | A x |^2
    scale.ht <- function (x, loc) {
        Ax <- G%*%x - dG*x
        Ax[loc] <- 0
        omitthese <- ( abs(scale(Ax,median(Ax),mad(Ax))) > 2 )
        Ax[omitthese] <- 0
        return( (-1)*sum(Ax)/sum(Ax^2) )
    }

    #  (d/dc) | A(x + c 1) - b |^2 = 2 1^T A^T ( A(x + c 1) - b )
    #    = 0  =>  c = - 1^T A^T (Ax-b) / 1^T A^T A 1
    shift.ht <- function (x, loc) {
        numerator <- G%*%x - dG*x + 1
        numerator[loc] <- 0
        numerator <- crossprod(G,numerator) - dG*numerator
        numerator[loc] <- 0
        denom <- sum( G[-loc,loc]^2 )
        return( (-1)*sum(numerator)/denom )
    }

    H.time <- system.time( sapply( 1:ncol(init.hts), function (k) H(init.hts[,k],loc=neighborhoods[[k]]) ) ) / ncol(init.hts)
    dH.time <- system.time( sapply( 1:ncol(init.hts), function (k) dH(init.hts[,k],loc=neighborhoods[[k]]) ) ) / ncol(init.hts)
    est.time <- ( maxit * (H.time[1] + dH.time[1]) / numcores * ncol(init.hts) )
    cat("Estimated time: ", est.time, " .\n")

    # maxit <- floor( maxtime / (H.time[1] + dH.time[1]) * numcores / ncol(init.hts) ) / 10  # turns out optim adds in a fair bit of overhead, hence the '/10'
    # maxit <- min( 1e4, maxit )

    optim.ht.list <- mclapply( seq_along(neighborhoods), function (loc.ind) {
                new.ht <- list(par=init.hts[,loc.ind])
                for (k in 1:10) {
                    aval <- if (k<nscale) { scale.ht(new.ht$par, neighborhoods[[loc.ind]]) } else { 1 }
                    bval <- if (k<nscale) { shift.ht(aval*new.ht$par, neighborhoods[[loc.ind]]) } else { 0 }
                    new.ht <- optim( par=pmax(0,aval*new.ht$par+bval), fn=H, gr=dH, loc=neighborhoods[[loc.ind]], 
                        method="L-BFGS-B", control=list( parscale=parscale, maxit=ceiling(maxit/10) ), lower=0, upper=Inf ) 
                    new.ht$par[neighborhoods[[loc.ind]]] <- 0
                }
                return(new.ht)
            }, mc.cores=numcores )

    convergences <- sapply(optim.ht.list,"[[","convergence")
    unconverged <- which(convergences != 0)

    # save( optim.ht.list, file=gsub(".tsv", "-optim.RData", outfile) )

    hts <- sapply( optim.ht.list, "[[", "par" )

}

colnames( hts ) <- locs
write.table( hts, file=outfile, row.names=FALSE )
cat("Writing output to ", outfile, " .\n")

## look at results
if (FALSE) {
        load( paste(subdir, "/", basename(layer.prefix),"_", basename(layer.file),"_nonmissing.RData",sep='') ) # provides nonmissing
        ph <- plot.ht.fn(layer.prefix,"annual_precip",nonmissing)

        layout(matrix(1:6,nrow=2))
        for (k in 1:ncol(hts)) {
            ph(hts[,k])
            if (k%%6 == 0 && is.null(locator(1))) break
        }
}

if (FALSE) {  ## DEBUGGING/EXPLORATION
    # check gradient
    loc.ind <- 10
    loc <- locs[loc.ind]
    ht <- init.hts[,loc.ind]

    eps <- 1e-5 * runif(length(ht))
    c( H(ht,loc=loc), H(ht+eps,loc=loc)-H(ht,loc=loc), sum(eps*dH(ht,loc=loc)) )

    bdry <- ( abs(dH(ht,loc=loc)) > 1e-4 )
    eps <- 1e-5 * ifelse( bdry, runif(length(ht)), 0 )
    c( H(ht,loc=loc), H(ht+eps,loc=loc)-H(ht,loc=loc), sum(eps*dH(ht,loc=loc)) )

    eps <- 1e-5 * ifelse( seq_along(ht) == sample(which(bdry),1), 1, 0 )
    c( H(ht,loc=loc), H(ht+eps,loc=loc)-H(ht,loc=loc), sum(eps*dH(ht,loc=loc)) )


    # look at convergence
    load( paste(subdir, "/", basename(layer.prefix),"_", basename(layer.file),"_nonmissing.RData",sep='') ) # provides nonmissing
    ph <- plot.ht.fn(layer.prefix,"annual_precip",nonmissing)


    loc.ind <- 10
    nreps <- 20
    parscale <- pmax(.01*mean(init.hts[,loc.ind]),init.hts[,loc.ind])
    oht.list <- vector( mode='list', length=nreps+1 )
    oht.list[[1]] <- list( par=init.hts[,loc.ind] )
    for (k in 2:(nreps+1)) { 
        oht.list[[k]] <- optim( par=oht.list[[k-1]]$par, fn=H, gr=dH, loc=neighborhoods[[loc.ind]], method="L-BFGS-B", control=list( parscale=parscale, maxit=100 ), lower=0, upper=Inf )  
    }

    oht.mat <- sapply(oht.list, "[[", "par")
    oht.dH.mat <- sapply( lapply(oht.list, "[[", "par"), dH, loc=neighborhoods[[loc.ind]])
    bdry <- ( abs(dH(oht.list[[1]]$par,loc=neighborhoods[[loc.ind]])) > quantile(abs(oht.dH.mat),.999) )
    plot.locs <- c( which(bdry), sample.int(nrow(oht.mat),100) )
    layout(matrix(1:4,nrow=2))
    matplot( t(oht.mat[plot.locs,]), type='l', lty=(1:2)[c(rep(1,sum(bdry)),rep(2,100))] )
    matplot( t((oht.mat-oht.mat[,1])[plot.locs,]), type='l', lty=(1:2)[c(rep(1,sum(bdry)),rep(2,100))] )
    matplot( t(oht.dH.mat[plot.locs,]), type='l', lty=(1:2)[c(rep(1,sum(bdry)),rep(2,100))] )
    matplot( abs(t(oht.dH.mat[plot.locs,])), type='l', lty=(1:2)[c(rep(1,sum(bdry)),rep(2,100))], log='y' )

    layout( matrix(1:6,nrow=2) )
    for (k in 2:length(oht.list)) { 
        ph( oht.list[[k]]$par )
        ph( oht.list[[k]]$par - oht.list[[k-1]]$par )
        ph( dH( oht.list[[k-1]]$par, loc=neighborhoods[[loc.ind]] ) )
        if (is.null(locator(1))) break
    }



    Gcheck <- function (x,loc.ind) { z <- as.vector(G%*%x-dG*x+1); z[neighborhoods[[loc.ind]]] <- 0; ph(z); invisible(z) }
    layout(t(1:2))
    aval <- scale.ht(init.hts[,loc.ind],neighborhoods[[loc.ind]]); 
    bval <- shift.ht(init.hts[,loc.ind],neighborhoods[[loc.ind]]); 
    plot( hts[,loc.ind], init.hts[,loc.ind] ); 
    abline(0,1); abline(0,1/aval,col='red')
    abline(-bval,1,col='green')
    Gcheck( init.hts[,loc.ind], loc.ind )
    Gcheck( aval*init.hts[,loc.ind], loc.ind )
    Gcheck( pmax(0,init.hts[,loc.ind]+bval), loc.ind )

    oht.list.1 <- vector( mode='list', length=nreps+1 )
    oht.list.1[[1]] <- list( par=init.hts[,loc.ind] )
    avec <- bvec <- numeric(nreps)
    for (k in 2:(nreps+1)) { 
        avec[k-1] <- if (TRUE) { scale.ht(oht.list.1[[k-1]]$par,neighborhoods[[loc.ind]]) } else { 1 }
        bvec[k-1] <- if (TRUE) { shift.ht(avec[k-1]*oht.list.1[[k-1]]$par,neighborhoods[[loc.ind]]) } else { 0 }
        oht.list.1[[k]] <- optim( par=pmax(0,avec[k-1]*oht.list.1[[k-1]]$par+bvec[k-1]), fn=H, gr=dH, loc=neighborhoods[[loc.ind]], method="L-BFGS-B", control=list( parscale=parscale, maxit=100 ), lower=0, upper=Inf )  
    }

    resids <- sapply(1:(nreps+1), function (k) { sum( ( hts[-neighborhoods[[loc.ind]],loc.ind] - oht.list.1[[k]]$par[-neighborhoods[[loc.ind]]] )^2 ) } )
    cbind(seq_along(resids),resids,c(0,avec),c(0,bvec))

    ##
    # THE GOOD STUFF:
    layout( matrix(1:6,nrow=2) )
    for (k in 2:length(oht.list.1)) { 
        ph( oht.list.1[[k]]$par - oht.list.1[[k-1]]$par, main=k )
        ph( dH( oht.list.1[[k-1]]$par, loc=neighborhoods[[loc.ind]] ) )
        plot( hts[-neighborhoods[[loc.ind]],loc.ind], oht.list.1[[k]]$par[-neighborhoods[[loc.ind]]] ); abline(0,1)
        points( hts[-neighborhoods[[loc.ind]],loc.ind], oht.list[[k]]$par[-neighborhoods[[loc.ind]]], col='green', pch=20 )
        if (is.null(locator(1))) break
    }

    oht.mat <- sapply(oht.list.1, "[[", "par")
    oht.dH.mat <- sapply( lapply(oht.list.1, "[[", "par"), dH, loc=neighborhoods[[loc.ind]])
    bdry <- ( abs(dH(oht.list.1[[1]]$par,loc=neighborhoods[[loc.ind]])) > quantile(abs(oht.dH.mat),.999) )
    plot.locs <- c( which(bdry), sample.int(nrow(oht.mat),100) )
    layout(matrix(1:4,nrow=2))
    matplot( t(oht.mat[plot.locs,]), type='l', lty=(1:2)[c(rep(1,sum(bdry)),rep(2,100))] )
    matplot( t((oht.mat-oht.mat[,1])[plot.locs,]), type='l', lty=(1:2)[c(rep(1,sum(bdry)),rep(2,100))] )
    matplot( t(oht.dH.mat[plot.locs,]), type='l', lty=(1:2)[c(rep(1,sum(bdry)),rep(2,100))] )
    matplot( abs(t(oht.dH.mat[plot.locs,])), type='l', lty=(1:2)[c(rep(1,sum(bdry)),rep(2,100))], log='y' )

    # overall difference
    oht.diff <- oht.list[[nreps]]$par-init.hts[,loc.ind]
    layout( matrix(1:6,nrow=2) )
    ph(oht.list[[nreps]]$par)
    ph(oht.diff)
    ph(log10(-oht.diff*(oht.diff<0)))
    ph(log10(oht.diff*(oht.diff>0)))
    ph(log10(abs(oht.diff/oht.list[[nreps]]$par)))

    layout( matrix(1:6,nrow=2) )
    ph( dH( oht.list[[length(oht.list)]]$par, loc=neighborhoods[[loc.ind]] ) )
    for (k in floor(seq(1,length(oht.list)-1,length.out=5))) {
        ph( oht.list[[length(oht.list)]]$par - oht.list[[k]]$par, main=k )
    }

    # look at inferred layer
    ph( (transfn(valfn(init.params[1 + (1:ngamma)]))) )
    ph( log(transfn(valfn(init.params[1 + (1:ngamma)]))) )
    ph( (transfn(valfn(init.params[1 + ngamma + (1:ndelta)]))) )
    ph( log(transfn(valfn(init.params[1 + ngamma + (1:ndelta)]))) )

}
