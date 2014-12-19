# cribbed from fit-logistic-model.R

layer.prefix <- "../geolayers/multigrid/256x/crm_"
subdir <- "256x"
layer.file <- "six-raster-list"
ht.file <- "test01/256x/six-raster-list-hitting-times-full.tsv"
param.file <- "test01/six-params.tsv"
outfile <- "test01/inferred-six-params.tsv"


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

layout( matrix( 1:16,nrow=4 ) )
plot.nearby( f=L, params=init.params, fac=.05 )

layout( matrix( 1:16,nrow=4 ) )
plot.nearby( f=L, params=init.params, fac=.5 )

fitfn <- function (init.params) {
    print(init.params)
    parscale <- c( abs(init.params[1]/10), rep(0.01,length(init.params)-1) )
    results <- optim( par=init.params, fn=L, gr=dL, control=list(parscale=parscale,fnscale=max(1,abs(L(init.params))/10)), method="BFGS" )
    # results <- optim( par=results$par, fn=L, gr=dL, control=list(parscale=parscale,fnscale=max(1,abs(L(results$par))/10)), method="BFGS" )
    results <- optim( par=results$par, fn=L, gr=dL, control=list(parscale=parscale,fnscale=max(1,abs(L(results$par))/10)), method="BFGS" )
    if (results$convergence != 0) {
        results <- optim( par=results$par, fn=L, gr=dL, control=list(parscale=parscale/10), method="BFGS" )
    }
}

jittered.params <- t( do.call( cbind, lapply( c(1e-3,1e-2,1e-1,1), function (sd) {  replicate( 2, init.params+rnorm(length(init.params),sd=sd) ) } ) ) )

fitted.params <- apply( jittered.params, 1, fitfn )
