#!/usr/bin/Rscript

if (!interactive()) {
    layer.prefix <- commandArgs(TRUE)[1]
    subdir <- commandArgs(TRUE)[2]
    layer.file <- commandArgs(TRUE)[3]
    ht.file <- commandArgs(TRUE)[4]
    param.file <- commandArgs(TRUE)[5]
} else {
    layer.prefix <- "../geolayers/multigrid/256x/crm_"
    subdir <- "256x"
    layer.file <- "dem-layer-list"
    ht.file <- "256x/dem-layer-list-hitting-times.tsv"
    param.file <- "params-dem-layer-list.tsv"
}
layer.names <- scan(layer.file,what="char") 

# require(raster)
# rasterOptions(tmpdir=".")

ht.id <- gsub( paste(basename(layer.prefix),basename(layer.file),sep=""), "", gsub("[.].*","", basename(ht.file) ) )
if (nchar(ht.id)>0) { ht.id <- paste("_",ht.id,sep='') }
outfile <- paste(subdir,"/",basename(layer.prefix),basename(layer.file),ht.id,"_fit-params.RData",sep='')

load(paste(subdir,"/",basename(layer.prefix),basename(layer.file),"-","setup.RData",sep=''))
hts <- as.matrix(read.table(ht.file,header=TRUE))

init.params.table <- read.table(param.file,header=TRUE)
init.params <- as.numeric( init.params.table[match(subdir,init.params.table[,1]),-1] )
names(init.params) <- colnames(init.params.table)[-1]

source("resistance-fns.R")

# Massage the numerics.
zeros <- unlist(neighborhoods) + rep((seq_along(neighborhoods)-1)*nrow(hts),sapply(neighborhoods,length))
scaling <- sqrt(nrow(G) * length(locs))
hts <- hts/scaling
hts[zeros] <- 0
usethese <- which( rowMeans(hts) < quantile(hts,.5) )  # indexes locations; note may overlap with zeros
nomitted <- ncol(hts)*length(usethese) + length( which(! row(hts)[zeros] %in% usethese) )
sc.one <- 1/scaling

# setup: evaluating L and dL will CHANGE THESE (once, for efficiency)
params <- init.params
G@x <- update.G(params)
dG <- rowSums(G)
GH <- G %*% hts - dG*hts
GH[zeros] <- 0

L <- function (params) {
    if (any(params != get("params",parent.env(environment()) ) ) ) { 
        assign("params", params, parent.env(environment()) )
        evalq( G@x <- update.G(params), parent.env(environment()) )
        evalq( dG <- rowSums(G), parent.env(environment()) )
        GH <- G %*% hts - dG*hts
        GH[zeros] <- 0
        assign("GH", GH, parent.env(environment()) )
    }
    ans <- ( sum( (GH+sc.one)[usethese,]^2 ) - (nomitted)*sc.one^2 )
    if (!is.finite(ans)) { browser() }
    return(ans)
}
dL <- function (params) {
    if (any(params != get("params", parent.env(environment()) ) ) ) { 
        assign("params", params, parent.env(environment()) )
        evalq( G@x <- update.G(params), parent.env(environment()) )
        evalq( dG <- rowSums(G), parent.env(environment()) )
        GH <- G %*% hts - dG*hts
        GH[zeros] <- 0
        assign("GH", GH, parent.env(environment()) )
    }
    bgrad <- ( 2 / params[1] )* sum( GH[usethese,] * (GH+sc.one)[usethese,] )
    ggrads <- sapply( 1:ncol(layers), function (kk) {
            2 * sum( (layers[,kk] * GH)[usethese,] * (GH+sc.one)[usethese,] )
        } )
    dgrads <- ggrads + sapply( 1:ncol(layers), function (kk) {
            GL <- G
            GL@x <- G@x * layers[Gjj,kk]
            dGL <- rowSums(GL)
            GLH <- GL %*% hts - dGL*hts
            GLH[zeros] <- 0
            return( 2 * sum( GLH[usethese,] * (GH+sc.one)[usethese,]  ) )
        } )
    ans <- ( c(bgrad, ggrads, dgrads) )
    if (any(!is.finite(ans))) { browser() }
    return(ans)
}
# environment(L) <- environment(dL) <- fun.env <- list2env( list(
#                 G=G,
#                 dG=dG,
#                 params=init.params,
#                 GH=GH), 
#         parent=environment() )

L(init.params)
dL(init.params)

parscale <- c( abs(init.params[1]/10), rep(0.1,length(init.params)-1) )
results <- optim( par=init.params, fn=L, gr=dL, control=list(parscale=parscale,fnscale=max(1,abs(L(init.params))/10)), method="BFGS" )
results <- optim( par=results$par, fn=L, gr=dL, control=list(parscale=parscale,fnscale=max(1,abs(L(results$par))/10)), method="BFGS" )
results <- optim( par=results$par, fn=L, gr=dL, control=list(parscale=parscale,fnscale=max(1,abs(L(results$par))/10)), method="BFGS" )
if (results$convergence != 0) {
    results <- optim( par=results$par, fn=L, gr=dL, control=list(parscale=parscale/10), method="BFGS" )
}

write( results$par, file=outfile )

if (FALSE) {

    # check gradient  (VERY STEEP ?!?!?!)
    dp <- 1e-8 * runif(length(init.params))
    L0 <- L(init.params)
    dL0 <- dL(init.params)
    L1 <- L(init.params+dp)
    dL1 <- dL(init.params+dp)
    c( L1-L0, sum(dp*dL0), sum(dp*dL1) )


    # check answer
    layout(t(seq_along(init.params)))
    for (k in seq_along(init.params)) {
        parvals <- seq( results$par[k]/1.02, results$par[k]*1.02, length.out=20 )
        Lvals <- sapply(parvals, function (x) L(ifelse(seq_along(init.params)==k,x,results$par)) )
        yrange <- range(Lvals,L(results$par))
        plot( parvals, Lvals, ylim=yrange, main=names(init.params)[k] )
        abline(v=results$par[k])
        abline(h=L(results$par))
    }

    G@x <- update.G(results$par)
    dG <- rowSums(G)
    GH <- G %*% hts - dG*hts
    GH[zeros] <- 0
    range(GH*scaling)

    # CHECK AGAINST RIGHT ANSWER
    ph <- plot.ht.fn(layer.prefix,"dem_30_m800_sq",nonmissing)
    params <- init.params
    G@x <- update.G(params)
    dG <- rowSums(G)
    hts <- hitting.analytic( neighborhoods, G-diag(rowSums(G)), numcores=getcores() )
    GH <- G %*% hts - dG*hts
    GH[zeros] <- 0

}
