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

# SWITCH TO LOGISTIC MODEL
transfn <- function (x) { 1/(1+exp(-x)) }

hts <- as.matrix(read.table(ht.file,header=TRUE))

init.params.table <- read.table(param.file,header=TRUE)
init.params <- as.numeric( init.params.table[match(subdir,init.params.table[,1]),-1] )
names(init.params) <- colnames(init.params.table)[-1]

source("resistance-fns.R")

# Massage the numerics.
zeros <- unlist(neighborhoods) + rep((seq_along(neighborhoods)-1)*nrow(hts),sapply(neighborhoods,length))
scaling <- 1 # sqrt(nrow(G) * length(locs))
sc.one <- 1/scaling
hts <- hts/scaling
hts[zeros] <- 0

# setup: evaluating L and dL will CHANGE THESE GLOBALLY (once, for efficiency)
update.aux <- function (params,env,check=TRUE) {
    if ( (!check) || any(params != get("params", env ) ) ) { 
        assign("params", params,  env )
        evalq( G@x <- update.G(params),  env )
        evalq( dG <- rowSums(G),  env )
        GH <- G %*% hts - dG*hts
        GH[zeros] <- 0
        assign("GH", GH, env )
    }
}

# weightings <- ifelse( rowMeans(hts) < quantile(hts,.5), dG, 0 )  # indexes locations; note may overlap with zeros
# weightings <- ifelse( 1:nrow(hts) %in% locs, 1, 0 )
weightings <- 1/rowMeans(hts,na.rm=TRUE)
nomitted <- sum( weightings[row(hts)[zeros]] )
update.aux(init.params,environment(),check=FALSE)

L <- function (params) {
    update.aux(params,parent.env(environment()))
    ans <- ( sum( weightings*rowSums((GH+sc.one)^2) ) - (nomitted)*sc.one^2 )
    if (!is.finite(ans)) { browser() }
    return(ans)
}
dL <- function (params) {
    update.aux(params,parent.env(environment()))
    gamma <- params[1+(1:ngamma)]
    bgrad <- ( 2 / params[1] )* sum( weightings * rowSums(GH * (GH+sc.one)) )
    ggrads <- sapply( 1:ncol(layers), function (kk) {
            2 * sum( weightings * rowSums( (layers[,kk] * (1-transfn(valfn(gamma))) * GH) * (GH+sc.one)) )
        } )
    dgrads <- ggrads + sapply( 1:ncol(layers), function (kk) {
                XXX EDITING HERE XXX
            GL <- G
            GL@x <- G@x * layers[Gjj,kk]
            dGL <- rowSums(GL)
            GLH <- GL %*% hts - dGL*hts
            GLH[zeros] <- 0
            return( 2 * sum( weightings * rowSums( GLH * (GH+sc.one) )  ) )
        } )
    ans <- ( c(bgrad, ggrads, dgrads) )
    if (any(!is.finite(ans))) { browser() }
    return(ans)
}

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

    gcheck <- function (params=jitter(init.params),eps=1e-8,dp=eps*dirn,dirn=runif(length(params))) {
        # check gradient  (VERY STEEP ?!?!?!)
        L0 <- L(params)
        dL0 <- dL(params)
        L1 <- L(params+dp)
        dL1 <- dL(params+dp)
        c( L0, L1-L0, sum(dp*dL0), sum(dp*dL1) )
    }

    # check answer
    layout(matrix(1:(2*length(init.params)),nrow=2,byrow=TRUE))
    for (fac in c(1.02,10)) {
        for (k in seq_along(init.params)) {
            parvals <- seq( results$par[k]/fac, results$par[k]*fac, length.out=20 )
            Lvals <- sapply(parvals, function (x) L(ifelse(seq_along(init.params)==k,x,results$par)) )
            yrange <- range(Lvals,L(results$par))
            plot( parvals, Lvals, ylim=yrange, main=names(init.params)[k] )
            abline(v=results$par[k])
            abline(h=L(results$par))
        }
    }

    update.aux(results$par)
    range(GH*scaling)

    # CHECK AGAINST RIGHT ANSWER
    ph <- plot.ht.fn(layer.prefix,"dem_30_m800_sq",nonmissing)
    params <- init.params
    G@x <- update.G(params)
    dG <- rowSums(G)
    hts <- hitting.analytic( neighborhoods, G-diag(rowSums(G)), numcores=getcores() )
    update.aux(params,environment(),check=FALSE)

    results <- optim( par=init.params*1.1, fn=L, gr=dL, control=list(parscale=parscale,fnscale=max(1,abs(L(init.params))/10)), method="BFGS" )

}
