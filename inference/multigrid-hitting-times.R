#!/usr/bin/Rscript

usage <- '
    Get hitting times with a list of landscape layers:
        Rscript multigrid-hitting-times.R (layer file) (parameter file) (initial guess) (output file) (dir.1) [dir.2] ...
      e.g.
        Rscript multigrid-hitting-times.R six-raster-list simple-init-params-six-raster-list.tsv multigrid-hts.tsv 

'

if (!interactive()) {
    if (length(commandArgs(TRUE))<8) { cat(usage) }
    layer.dir <- commandArgs(TRUE)[1]
    layer.prefix <- commandArgs(TRUE)[2]
    layer.file <- commandArgs(TRUE)[3]
    param.file <- commandArgs(TRUE)[4] 
    prev.ht <- commandArgs(TRUE)[5]
    outfile <- commandArgs(TRUE)[6]
    dirnames <- commandArgs(TRUE)[-(1:6)]
} else {
    layer.dir <- "../geolayers/multigrid/"
    layer.prefix <- "crm_"
    layer.file <- "six-raster-list"
    param.file <- "simple-init-params-six-raster-list.tsv"
    prev.ht <- "512x/six-raster-list-hitting-times.tsv"
    outfile <- "temp-hitting.tsv"
    dirnames <- c("512x","256x","128x")
}
names(dirnames) <- dirnames

ag.fact <- 2

source("resistance-fns.R")
require(raster)

require(parallel)
numcores<-getcores()

layer.names <- scan(layer.file,what="char") 

res.envs <- lapply( dirnames, function (subdir) { 
        env <- new.env(parent=.GlobalEnv)
        assign("subdir",subdir,env)
        #   provides "G"    "Gjj"    "update.G" "ndelta"   "ngamma"   "transfn"  "valfn"    "layers"
        load( paste(subdir,"/",basename(layer.prefix),"_",basename(layer.file),"_","G.RData",sep=''), envir=env ) 
        # provides 'neighborhoods'
        load( paste( subdir, "/", basename(layer.prefix), basename(layer.file), "_neighborhoods.RData", sep='' ), envir=env ) 
        # provides 'locs'
        load(paste(subdir,"/",basename(layer.prefix),basename(layer.file),"_tortlocs.RData",sep=''), envir=env) 
        # REMOVE MISSING INDIV
        na.indiv <- which( is.na( get("locs",env) ) )
        assign("locs", get("locs",env)[-na.indiv], env) 
        assign("neighborhoods", lapply(get("neighborhoods",env)[-na.indiv],function (x) { x[!is.na(x)] }), env)
        assign("dG", rowSums(get("G",env)), env) 
        assign("layer", raster(paste(file.path(layer.dir,subdir,layer.prefix),"annual_precip",sep='')), env)
        # provides nonmissing
        load( paste(subdir, "/", basename(layer.prefix),"_", basename(layer.file),"_nonmissing.RData",sep=''), envir=env ) 
        return(env)
    } )


##
#  Parameters
#   in time t, 1D RW with rate r does t*r jumps, displacement has variance t*r
#   so time to move N grid sites away is sqrt(N)/r
#   so if hitting times are of order T, want r of order sqrt(N)/T
init.param.table <- read.table( param.file, header=TRUE )
for (k in 1:nrow(init.param.table)) {
    subdir <- init.param.table[k,1]
    if (subdir %in% names(res.envs)) {
        assign("G@x", update.G(init.param.table[k,-1]),res.envs[[subdir]])
    }
}


##
# Initial hitting times
orig.init.hts <- as.matrix( read.table(prev.ht,header=TRUE) )
init.hts <- sapply( 1:ncol(orig.init.hts), function (k) {
            x <- orig.init.hts[,k]
            x[get("neighborhoods",res.envs[[1]])[[k]]] <- 0  # these are from FIRST resolution
            return(x) } )


# objective function
H <- function (ht,loc.ind) {
    # ( (G-D)ht + 1 )^T S ( (G-D)ht + 1)
    # where S = I except S[loc,loc]=0
    loc <- neighborhoods[[loc.ind]]
    ht[loc] <- 0
    z <- G%*%ht - dG*ht + 1
    z[loc] <- 0
    return( ( sum( z^2 ) )/length(z) )
}
dH <- function (ht,loc.ind) {
    # 2 (G-D)^T S ( (G-D) ht + 1 )
    loc <- neighborhoods[[loc.ind]]
    ht[loc] <- 0
    z <- G%*%ht - dG*ht + 1
    z[loc] <- 0
    z <- (crossprod(G,z) - dG*z)
    z[loc] <- 0
    return( 2 * as.vector(z) / length(z) )
}

loc.ind <- 10

stepupdn <- function (init.hts,oldpos,newpos) {
    return( if (newpos>oldpos) {
        upsample( init.hts, ag.fact, 
            get("layer",res.envs[[oldpos]]), get("nonmissing",res.envs[[oldpos]]), 
            get("layer",res.envs[[newpos]]), get("nonmissing",res.envs[[newpos]]) )
    } else if (newpos<oldpos) {
        downsample( init.hts, ag.fact, 
            get("layer",res.envs[[newpos]]), get("nonmissing",res.envs[[newpos]]), 
            get("layer",res.envs[[oldpos]]), get("nonmissing",res.envs[[oldpos]]) )
    } else {
        init.hts
    } )
}

step <- function (init.hts,oldpos,newpos,maxit) {
    ihts <- stepupdn(init.hts,oldpos,newpos)
    environment(H) <- environment(dH) <- res.envs[[newpos]]

    parscale <- (mean(ihts)+ihts)/2

    # optim.ht.list <- mclapply( seq_along(neighborhoods), function (loc.ind) {
            optim( par=ihts, fn=H, gr=dH, loc=loc.ind, 
                method="L-BFGS-B", control=list( parscale=parscale, maxit=maxit ), lower=0, upper=Inf ) 
    #   }, mc.cores=numcores )
    # return(sapply( optim.ht.list, "[[", "par" ))
}

# for plotting
phs <- lapply( res.envs, function (env) { with(env, plot.ht.fn(file.path(layer.dir,subdir,layer.prefix),"annual_precip",nonmissing)) } )

stepres <- rep(c(1,2),20)
maxits <- c(10,100)[stepres[-1]]
steps <- vector( length(stepres), mode="list" )
steps[[1]] <- list( par=init.hts[,loc.ind] )
for (k in seq_along(stepres)[-1]) {
    cat(k,"\n")
    steps[[k]] <- step( steps[[k-1]]$par, stepres[k-1], stepres[k], maxits[k-1] )
}

# differences
layout(matrix(1:6,nrow=2))
for (k in seq_along(stepres)[-1]) {
    prev <- stepupdn( steps[[k-1]]$par, stepres[k-1], stepres[k] )
    phs[[stepres[k]]](steps[[k]]$par,main=paste(k,":",dirnames[stepres[k]]))
    phs[[stepres[k]]](steps[[k]]$par-prev,main=paste(k,":",dirnames[stepres[k]]))
    environment(dH) <- res.envs[[stepres[k]]]
    phs[[stepres[k]]](dH(steps[[k]]$par,loc.ind),main=paste(k,":",dirnames[stepres[k]]))
    environment(dH) <- .GlobalEnv
    if (is.null(locator(1))) break
}

hts.mat <- sapply( steps[stepres==2], "[[", "par" )
environment(dH) <- res.envs[[2]]
dH.mat <- sapply( lapply( steps[stepres==2], "[[", "par" ), dH, loc=loc.ind)
environment(dH) <- .GlobalEnv
bdry <- ( abs(dH.mat[,1]) > quantile(abs(dH.mat),.995) )
plot.locs <- c( which(bdry), sample.int(nrow(hts.mat),60) )
layout(matrix(1:4,nrow=2))
matplot( t(hts.mat[plot.locs,]), type='l', lty=(1:2)[c(rep(1,sum(bdry)),rep(2,100))] )
matplot( t((hts.mat-hts.mat[,1])[plot.locs,]), type='l', lty=(1:2)[c(rep(1,sum(bdry)),rep(2,100))] )
matplot( t(dH.mat[plot.locs,]), type='l', lty=(1:2)[c(rep(1,sum(bdry)),rep(2,100))] )
matplot( abs(t(dH.mat[plot.locs,])), type='l', lty=(1:2)[c(rep(1,sum(bdry)),rep(2,100))], log='y' )
