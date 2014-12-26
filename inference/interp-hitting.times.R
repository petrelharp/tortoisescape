#!/usr/bin/Rscript

usage <- "
     Interpolate provided hitting times based on parameters given
         Rscript initial-hitting-times.R (layer prefix) (subdir) (layer file) (parameter file) (hitting times) (alpha) (output file)
       e.g.
         Rscript initial-hitting-times.R ../geolayers/multigrid/256x/crm_ 256x six-raster-list test_six_layers/six-params.tsv test_six_layers/256x/six-raster-list-sim-0_00-hts.tsv 1.0 test_six_layers/256x/six-raster-list-interp-hts.tsv 

 "

argvec <- if (!interactive()) { commandArgs(TRUE) } else { scan(what='char') }
if (length(argvec)<6) { stop(usage) }

layer.prefix <- argvec[1]
subdir <- argvec[2]
layer.file <- argvec[3]
param.file <- argvec[4]
ht.file <- argvec[5]
alpha <- as.numeric( argvec[6] )
output.file <- argvec[7]

cat("interp-hitting-times.R:\n")
invisible( lapply( c("layer.prefix","subdir","layer.file","param.file","ht.file","alpha","output.file"), function (x) { cat("  ", x, " : ", get(x), "\n") } ) )
cat("\n")

layer.names <- scan(layer.file,what="char") 

load(paste(subdir,"/",basename(layer.prefix),basename(layer.file),"-","setup.RData",sep=''))

source("resistance-fns.R")

require(parallel)
numcores<-getcores()

# hitting times
if (FALSE) {
    orig.obs.ht <- read.sub.hts( "test_six_layers/256x/six-raster-list-sim-0_00-hts.tsv", locs )
    obs.ht <- orig.obs.ht * exp( rnorm(length(orig.obs.ht)) * 1e-8 )
}

obs.ht.df <- read.table(ht.file,header=TRUE,stringsAsFactors=FALSE)
obs.ht <- matrix( NA, nrow=length(locs), ncol=length(locs) )
obs.ht[ cbind( obs.ht.df$row, obs.ht.df$col ) ] <- obs.ht.df$DISTANCE
# indices of these in the big matrix of hitting times is ht[locs,]

##
# initial parameters?
init.params.table <- read.table(param.file,header=TRUE)
init.params <- as.numeric( init.params.table[match(subdir,init.params.table[,1]),-1] )
names(init.params) <- colnames(init.params.table)[-1]

G@x <- update.G(init.params)
dG <- rowSums(G)

if (method=="analytic") {

    hts <- interp.hitting( neighborhoods, G-diag(rowSums(G)), obs.ht, obs.locs=locs, alpha=alpha, numcores=numcores )

    hts <- interp.hitting( neighborhoods[1:2], G-diag(rowSums(G)), obs.ht[,1:2], obs.locs=locs, alpha=alpha, numcores=numcores )

    interp.tradeoff( hts, neighborhoods[1:2], G, dG, obs.ht[,1:2], obs.locs=locs, numcores=numcores )

} else if (method=="numeric") {

    HdH <- interp.ht.setup()

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

###
# Conjugate gradient

dG <- rowSums(G)
cG <- colSums(G)
# objective function
H <- function (par,obs.ht,loc,locs,g.match=1) {
    a <- par[1]
    hts <- par[-1]
    hts[loc] <- 0
    z <- G%*%hts - dG*hts + 1
    z[loc] <- 0
    return( ( sum( z^2 ) + g.match * sum( (hts[locs] - (obs.ht-a) )^2 ) )/length(z) )
}
dH <- function (par,obs.ht,loc,locs,g.match=1) {
    # cG - G[loc,] is, except at [loc], 1^T ((G-diag(dG))[-loc,])
    a <- par[1]
    hts <- par[-1]
    z <- G%*%hts - dG*hts   # NOTE: Should be +1 here?
    z[loc] <- 0
    z <- (G%*%z - dG*z) + (cG-G[loc,])   # and without this last bit?
    z[locs] <- z[locs] + g.match*(hts[locs]-(obs.ht-a))
    z[loc] <- 0
    return( c( 2 * g.match * sum( hts[locs] - (obs.ht-a) ) / length(z) , 2 * as.vector(z) / length(z) ) )
}

# parscale <- rep( nrow(G) / exp( mean( log(dG), trim=.1, na.rm=TRUE ) ), nrow(G) )
parscale <- c( min(pimat), rep( mean(pimat), nrow(G) ) )
init.hts <- matrix(parscale,nrow=nrow(G)+1,ncol=length(locs))
init.hts[cbind(locs,seq_along(locs))] <- 0

# H( init.hts[,1], obs.ht=pimat[,1], loc=locs[1], locs=locs )
# dH( init.hts[,1], obs.ht=pimat[,1], loc=locs[1], locs=locs )

# k <- 28 
# test.ht <- optim( par=init.hts[,k], fn=H, gr=dH, obs.ht=pimat[,k], loc=locs[k], locs=locs, g.match=1/200,
#         method="L-BFGS-B", control=list( parscale=parscale, maxit=1000 ), lower=0, upper=Inf ) 
# ph(test.ht$par[-1]); with( environment(ph), points(tort.coords.rasterGCS[k]) )

optim.ht.list <- mclapply( seq_along(locs), function (k) {
            optim( par=init.hts[,k], fn=H, gr=dH, obs.ht=pimat[,k], loc=locs[k], locs=locs, 
                method="L-BFGS-B", control=list( parscale=parscale, maxit=1000 ), lower=0, upper=Inf ) 
        }, mc.cores=numcores )

convergences <- sapply(optim.ht.list,"[[","convergence")
unconverged <- which(convergences != 0)

for (k in 1:3) {
    if (length(unconverged)<length(locs)) {
        optim.ht.list[unconverged] <- mclapply( unconverged, function (k) {
                    newstart <- sample( setdiff(seq_along(locs),unconverged), 1 )
                    optim( par=optim.ht.list[[newstart]]$par, fn=H, gr=dH, obs.ht=pimat[,k], loc=locs[k], locs=locs, 
                        method="L-BFGS-B", control=list( parscale=parscale, maxit=1000 ), lower=0, upper=Inf ) 
                }, mc.cores=numcores )
        convergences <- sapply(optim.ht.list,"[[","convergence")
        unconverged <- which(convergences != 0)
    }
}

optim.hts <- sapply(optim.ht.list,"[[","par")

if (any(convergences!=0)) { warning("Some did not converge") }

save( layer.prefix, layer.names, subdir, optim.hts, convergences, init.params, ratescale, gridwidth, file=paste(subdir,"/",basename(layer.prefix),basename(layer.file),"-init-hts.RData",sep='') )


###
# testing
if (FALSE) {


    load(paste(subdir,"/",basename(layer.file),"-",basename(layer.prefix),"setup.RData",sep=''))
    load("../tort.coords.rasterGCS.Robj")
    load(paste(subdir,"/",basename(layer.prefix),"alllocs.RData",sep='')) # all.locs.dists
    load("torts-info.RData")
    ph <- plot.ht.fn(layer.prefix,"annual_precip",nonmissing)

    optim.hts[cbind(1+locs,seq_along(locs))] <- NA

    for (k in 1:ncol(optim.hts)) {
        layout(matrix(1:4,nrow=2))
        ph( optim.hts[-1,k], main=optim.hts[1,k] )
        points(tort.coords.rasterGCS[k])
        plot( all.locs.dists[,k], optim.hts[-1,k], pch=20, cex=.5 )
        points( all.locs.dists[locs,k], optim.hts[1+locs,k], col='red' )
        plot( tort.dists[k,], optim.hts[(1+locs),k], pch=20, cex=.5 )
        plot( pimat[,k]-optim.hts[1,k], optim.hts[(1+locs),k], pch=20, cex=.5 ); abline(0,1)
        if (is.null(locator(1))) { break }
    }


    ###
    # analytic

    tmp.pimat <- pimat-mean(diag(pimat))

    solve.hts <- interp.hitting( locs, fullG, tmp.pimat )
    solve.hts[cbind(locs,seq_along(locs))] <- 0

    plot( as.vector(tmp.pimat), as.vector(solve.hts[locs,]), col=1+(row(pimat)==col(pimat)) ); abline(0,1)

    for (k in seq_along(locs)[1:length(locs)]) {
        plot.ht( pmax(solve.hts[,k],0), hitting.layer, nonmissing )
        text( tort.coords.rasterGCS, labels=1:180 )
        points( tort.coords.rasterGCS[k+if(k>56){1}else{0}], pch="*", cex=4, col='red' )
        if (is.null(locator(1))) { break }
    }
}

