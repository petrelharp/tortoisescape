#!/usr/bin/Rscript

usage <- "
Fit the parameters in a model where the transformation function is the logistic function x -> 1/(1+exp(-x)).
Usage:
    Rscript fit-logistic-model.R (layer prefix) (subdir) (layer file) (hitting times tsv) (initial params file) [output file]
e.g.
    Rscript fit-logistic-model.R ../geolayers/multigrid/256x/crm_ 256x dem-layer-list 256x/dem-layer-list-hitting-times.tsv params-dem-layer-list.tsv
"

if (!interactive()) {
    layer.prefix <- commandArgs(TRUE)[1]
    subdir <- commandArgs(TRUE)[2]
    layer.file <- commandArgs(TRUE)[3]
    ht.file <- commandArgs(TRUE)[4]
    param.file <- commandArgs(TRUE)[5]
    outfile <- if (length(commandArgs(TRUE))>5) { commandArgs(TRUE)[6] } else { NULL }
} else {
    layer.prefix <- "../geolayers/multigrid/256x/crm_"
    subdir <- "256x"
    layer.file <- "dem-layer-list"
    ht.file <- "256x/dem-layer-list-hitting-times.tsv"
    param.file <- "params-dem-layer-list.tsv"
    outfile <- NULL
}

layer.names <- scan(layer.file,what="char") 

# require(raster)
# rasterOptions(tmpdir=".")

ht.id <- gsub( paste(basename(layer.prefix),basename(layer.file),sep=""), "", gsub("[.].*","", basename(ht.file) ) )
if (nchar(ht.id)>0) { ht.id <- paste("_",ht.id,sep='') }
if (is.null(outfile)) { outfile <- paste(subdir,"/",basename(layer.prefix),basename(layer.file),ht.id,"_fit-params.RData",sep='') }

load(paste(subdir,"/",basename(layer.prefix),basename(layer.file),"-","setup.RData",sep=''))

# MAKE SURE USING LOGISTIC FUNCTION
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

LdL <- params.logistic.setup()
L <- LdL$L
dL <- LdL$dL

parscale <- c( abs(init.params[1]/10), rep(0.1,length(init.params)-1) )
results <- optim( par=init.params, fn=L, gr=dL, control=list(parscale=parscale,fnscale=max(1,abs(L(init.params))/10)), method="BFGS" )
results <- optim( par=results$par, fn=L, gr=dL, control=list(parscale=parscale,fnscale=max(1,abs(L(results$par))/10)), method="BFGS" )
if (results$convergence != 0) {
    results <- optim( par=results$par, fn=L, gr=dL, control=list(parscale=parscale/10), method="BFGS" )
}

write( results$par, file=outfile )

if (FALSE) {

    gcheck <- function (params=jitter(init.params),eps=1e-8,dp=eps*dirn,dirn=runif(length(params))) {
        # check gradient: want L1-L0 == dL0 if eps is small enough that dL0 == dL1
        L0 <- L(params)
        dL0 <- dL(params)
        L1 <- L(params+dp)
        dL1 <- dL(params+dp)
        c( L0=L0, diffL=L1-L0, dL0=sum(dp*dL0), dL1=sum(dp*dL1) )
    }

    for (k in seq_along(init.params)) {
        cat("Checking ", names(init.params)[k], " .\n")
        print( gcheck(dirn=ifelse(seq_along(init.params)==k,1,0)) )
    }

    # check answer: marginal plots to +/- 0.1
    layout(matrix(1:(2*ceiling(length(init.params)/2)),nrow=2,byrow=TRUE))
    plot.nearby(L,results$par,fac=0.1)

    update.aux(results$par)
    range(GH*scaling)

    # CHECK AGAINST RIGHT ANSWER
    require(raster)
    ph <- plot.ht.fn(layer.prefix,"dem_30_m800_sq",nonmissing)
    params <- init.params
    G@x <- update.G(params)
    dG <- rowSums(G)
    hts <- hitting.analytic( neighborhoods, G-diag(rowSums(G)), numcores=getcores() )
    update.aux(params,check=FALSE)

    results <- optim( par=init.params*1.1, fn=L, gr=dL, control=list(parscale=parscale,fnscale=max(1,abs(L(init.params))/10)), method="BFGS" )

}
