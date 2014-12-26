# cribbed from fit-logistic-model.R

layer.prefix <- "../geolayers/multigrid/256x/crm_"
subdir <- "256x"
layer.file <- "six-raster-list"
ht.file <- "test_six_layers/256x/six-raster-list-hitting-times-full.tsv"
param.file <- "test_six_layers/six-params.tsv"
outfile <- "test_six_layers/inferred-six-params.tsv"


layer.names <- scan(layer.file,what="char") 

# require(raster)
# rasterOptions(tmpdir=".")

ht.id <- gsub( paste(basename(layer.prefix),basename(layer.file),sep=""), "", gsub("[.].*","", basename(ht.file) ) )
if (nchar(ht.id)>0) { ht.id <- paste("_",ht.id,sep='') }
if (is.null(outfile)) { outfile <- paste(subdir,"/",basename(layer.prefix),basename(layer.file),ht.id,"_fit-params.RData",sep='') }

load(paste(subdir,"/",basename(layer.prefix),basename(layer.file),"-","setup.RData",sep=''))

# MAKE SURE USING LOGISTIC FUNCTION
transfn <- function (x) { 1/(1+exp(-x)) }

# load hitting times
hts <- orig.hts <- as.matrix(read.table(ht.file,header=TRUE))
# add some noise to them
hts <- hts * exp( rnorm(length(hts))/100 )

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

parscale <- c( .01, rep(1,length(init.params)-1) )

# parscale <- find.parscale( L, init.params, parscale=rep(1,length(init.params)) )
# parscale <- c(1,rep(1,ngamma),rep(0.0001,ndelta))
check.parscale(L,init.params,parscale)

fitfn <- function (init.params) {  # takes like 5min if not close to true params
    print(init.params)
    results <- optim( par=init.params, fn=L, gr=dL, control=list(parscale=parscale,fnscale=max(1,abs(L(init.params))/10)), method="BFGS" )
    if (results$convergence != 0) {
        results <- optim( par=results$par, fn=L, gr=dL, control=list(parscale=parscale/10), method="BFGS" )
    }
    return(results)
}

results <- fitfn( init.params )

layout( matrix( 1:16,nrow=4 ) )
plot.nearby( f=L, params=results$par, fac=.5 )

jittered.params <- t( do.call( cbind, lapply( c(1e-3,1e-2,1e-1,1), function (sd) {  replicate( 2, init.params+rnorm(length(init.params),sd=sd) ) } ) ) )

fitted.params <- apply( jittered.params, 1, fitfn )


