#!/usr/bin/Rscript

require(raster)
rasterOptions(tmpdir=".")

if (!interactive()) {
    layer.prefix <- commandArgs(TRUE)[1]
    subdir <- commandArgs(TRUE)[2]
    layer.file <- commandArgs(TRUE)[3]
    ht.file <- if (length(commandArgs(TRUE))>3) { commandArgs(TRUE)[4] } else { paste(subdir,"/",basename(layer.prefix),basename(layer.file),"-init-hts.RData",sep='') }
} else {
    layer.prefix <- c("../geolayers/TIFF/500x/500x_")
    subdir <- "500x"
    layer.file <- "six-raster-list"
    ht.file <- paste(subdir,"/",basename(layer.prefix),basename(layer.file),"-init-hts.RData",sep='')
    # layer.names <- c("imperv_30", "agp_250", "m2_ann_precip", "avg_rough_30", "dem_30", "bdrock_ss2_st")
}
layer.names <- scan(layer.file,what="char") 

ht.id <- gsub( paste(basename(layer.prefix),basename(layer.file),sep=""), "", gsub("-init-hts.RData","", basename(ht.file) ) )
if (nchar(ht.id)>0) { ht.id <- paste("_",ht.id) }
outfile <- paste(subdir,"/",basename(layer.prefix),basename(layer.file),ht.id,"-fit-params.RData",sep='')

load(paste(subdir,"/",basename(layer.prefix),basename(layer.file),"-","setup.RData",sep=''))
load(ht.file)

source("resistance-fns.R")

hts <- optim.hts[-1,]
ht.shifts <- optim.hts[1,]
rm(optim.hts)

# Massage the numerics.
zeros <- which( row(hts) == locs[col(hts)] )
scaling <- sqrt(nrow(G) * length(locs))
hts <- hts/scaling
hts[zeros] <- 0
sc.one <- 1/scaling

#
dG <- rowSums(G)
GH <- G %*% hts - dG*hts
GH[zeros] <- 0

init.beta <- init.params[1]

L <- function (params) {
    if (any(params != get("params",parent.env(environment()) ) ) ) { 
        assign("params", params, parent.env(environment()) )
        evalq( G@x <- update.G(c(init.beta,params)), parent.env(environment()) )
        evalq( dG <- rowSums(G), parent.env(environment()) )
        GH <- G %*% hts - dG*hts
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
        GH <- G %*% hts - dG*hts
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
            GLH <- GL %*% hts - dGL*hts
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

write( c(init.beta,results$par), file=outfile )

if (FALSE) {

    # check gradient  (VERY STEEP ?!?!?!)
    dp <- 1e-8 * runif(length(init.params)-1)
    L0 <- L(init.params[-1])
    dL0 <- dL(init.params[-1])
    L1 <- L(init.params[-1]+dp)
    c( L1-L0, sum(dp*dL0) )



    # version with beta free to vary
L <- function (params) {
    if (any(params != get("params",parent.env(environment()) ) ) ) { 
        evalq(params <- params, parent.env(environment()) )
        evalq( G@x <- update.G(params), parent.env(environment()) )
        evalq( dG <- rowSums(G), parent.env(environment()) )
        GH <- G %*% hts - dG*hts
        GH[zeros] <- 0
        evalq(GH <- GH, parent.env(environment()) )
    }
    ans <- ( sum( (GH+sc.one)^2 ) - length(zeros)*sc.one^2 )
    if (!is.finite(ans)) { browser() }
    return(ans)
}
dL <- function (params) {
    if (any(params != get("params", parent.env(environment()) ) ) ) { 
        evalq(params <- params, parent.env(environment()) )
        evalq( G@x <- update.G(params), parent.env(environment()) )
        evalq( dG <- rowSums(G), parent.env(environment()) )
        GH <- G %*% hts - dG*hts
        GH[zeros] <- 0
        evalq(GH <- GH, parent.env(environment()) )
    }
    bgrad <- 2 * sum( GH * (GH+sc.one) ) / params[1]
    ggrads <- sapply( 1:ncol(layers), function (kk) {
            2 * sum( layers[,kk] * GH * (GH+sc.one) )
        } )
    dgrads <- ggrads + sapply( 1:ncol(layers), function (kk) {
            GL <- G
            GL@x <- G@x * layers[Gjj,kk]
            dGL <- rowSums(GL)
            GLH <- GL %*% hts - dGL*hts
            GLH[zeros] <- 0
            return( 2 * sum( GLH * (GH+sc.one)  ) )
        } )
    ans <- ( c(bgrad, ggrads, dgrads) )
    if (any(!is.finite(ans))) { browser() }
    return(ans)
}
environment(L) <- environment(dL) <- fun.env <- list2env( list(
                G=G,
                dG=dG,
                params=init.params,
                GH=GH), 
        parent=environment() )
parscale <- c( init.params[1]/10, rep(0.1,length(init.params)-1) )
results <- optim( par=init.params, fn=L, gr=dL, control=list(parscale=parscale, trace=5), method="L-BFGS-B", lower=c(1e-6,rep(-Inf,length(init.params)-1)) )

}
