#!/usr/bin/Rscript

usage <- "
     Interpolate provided hitting times based on parameters given
         Rscript initial-hitting-times.R (layer prefix) (subdir) (layer file) (parameter file) (observed hitting times file) (alpha) (method) (output file) [initial hitting times file] 
       e.g.
         Rscript initial-hitting-times.R ../geolayers/multigrid/256x/crm_ 256x six-raster-list test_six_layers/six-params.tsv test_six_layers/256x/six-raster-list-sim-0_00-hts.tsv 1.0 analytic test_six_layers/256x/six-raster-list-interp-hts.tsv 

 "

argvec <- if (!interactive()) { commandArgs(TRUE) } else { scan(what='char') }
if (length(argvec)<8) { stop(usage) }

layer.prefix <- argvec[1]
subdir <- argvec[2]
layer.file <- argvec[3]
param.file <- argvec[4]
ht.file <- argvec[5]
alpha <- as.numeric( argvec[6] )
method <- argvec[7]
output.file <- argvec[8]
init.ht.file <- if (length(argvec) > 8) { argvec[9] } else { NULL }
if ( (method!="analytic") && is.null(init.ht.file) ) { stop("With method 'numeric' need a file of initial hitting times.") }

maxit <- 100

cat("interp-hitting-times.R:\n")
invisible( lapply( c("layer.prefix","subdir","layer.file","param.file","ht.file","alpha","method","output.file","init.ht.file"), function (x) { cat("  ", x, " : ", get(x), "\n") } ) )
cat("\n")

layer.names <- scan(layer.file,what="char") 

load(paste(subdir,"/",basename(layer.prefix),basename(layer.file),"-","setup.RData",sep=''))

source("resistance-fns.R")

require(parallel)
numcores<-getcores()

if (FALSE) { # testing
    orig.obs.ht <- read.sub.hts( "test_six_layers/256x/six-raster-list-sim-0_00-hts.tsv", locs )
    obs.ht <- orig.obs.ht * exp( rnorm(length(orig.obs.ht)) * 1e-8 )
}

# hitting times
obs.ht <- read.sub.hts( ht.file, locs )
# indices of these in the big matrix of hitting times is ht[locs,]

##
# initial parameters?
init.params.table <- read.table(param.file,header=TRUE)
init.params <- as.numeric( init.params.table[match(subdir,init.params.table[,1]),-1] )
names(init.params) <- colnames(init.params.table)[-1]

G@x <- update.G(init.params)
dG <- rowSums(G)

if (method=="analytic") {

    # if alpha is small, finds hitting times for G; if alpha is large, finds values that are close to obs.ht;
    # as alpha increases may get large negative values
    # in places where migration rates are very small
    # ... how to pick alpha?

    hts <- interp.hitting( neighborhoods, G, obs.ht, obs.locs=locs, alpha=alpha, numcores=numcores )


    if (FALSE) { # testing
        hts.list <- lapply( seq(0.0,0.2,length.out=5), function (alpha) { 
                interp.hitting( neighborhoods[1:2], G, obs.ht[,1:2], obs.locs=locs, alpha=alpha, numcores=numcores )
            } )
        sapply( hts.list, interp.tradeoff, neighborhoods[1:2], G, dG, obs.ht[,1:2], obs.locs=locs, numcores=numcores )
    }

} else if (method=="numeric") {

    init.hts <- read.full.hts( init.ht.file, locs )

    HdH <- interp.ht.setup()
    H <- HdH$H
    dH <- HdH$dH

    parscale <- rep(mean(init.hts[locs,]),nrow(init.hts))

    optim.ht.list <- mclapply( seq_along(neighborhoods), function (loc.ind) {
                optim( par=init.hts[,loc.ind], fn=H, gr=dH, 
                        locs=locs, zeros=neighborhoods[[loc.ind]], alpha=alpha, obs.ht=obs.ht[,loc.ind],
                        method="L-BFGS-B", control=list( parscale=parscale, maxit=maxit, fnscale=sum(init.hts[,1]^2) ), lower=0, upper=Inf ) 
            }, mc.cores=numcores )

    convergences <- sapply(optim.ht.list,"[[","convergence")
    unconverged <- which(convergences != 0)

    optim.ht.list[unconverged] <- mclapply( seq_along(neighborhoods)[unconverged], function (loc.ind) {
                optim( par=optim.ht.list[[loc.ind]]$par, fn=H, gr=dH, 
                        locs=locs, zeros=neighborhoods[[loc.ind]], alpha=alpha, obs.ht=obs.ht[,loc.ind],
                        method="L-BFGS-B", control=list( parscale=parscale, maxit=maxit, fnscale=sum(init.hts[,1]^2) ), lower=0, upper=Inf ) 
            }, mc.cores=numcores )

    hts <- sapply( optim.ht.list, "[[", "par" )

}

write.full.hts( hts, locs, file=outfile )

##################


###
# testing
if (FALSE) {


    load(paste(subdir,"/",basename(layer.file),"-",basename(layer.prefix),"setup.RData",sep=''))
    load("../tort.coords.rasterGCS.Robj")
    load(paste(subdir,"/",basename(layer.prefix),"alllocs.RData",sep='')) # all.locs.dists
    load("torts-info.RData")
    ph <- plot.ht.fn(layer.prefix,nonmissing,"dem_30")

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

