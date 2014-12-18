#!/usr/bin/Rscript

require(raster)
source("resistance-fns.R")
require(parallel); numcores <- getcores()

biglayer <- raster(nrows=2^10,ncols=2^10,xmn=0,ymn=0,xmx=100,ymx=100)
values(biglayer) <- rnorm(length(biglayer))
nlayers <- 2
fact <- 2^5
layer.list <- lapply( 1:nlayers, function (k) {
        values(biglayer) <- rnorm(length(biglayer))
        aggregate( biglayer, fact=2^5, fun=mean, na.rm=FALSE )
    } )
navals <- sample(length(values(layer.list[[1]])),10)
for (k in 1:nlayers) { values(layer.list[[k]])[navals] <- NA }
layer <- layer.list[[1]]

nlocs <- 40
locs <- SpatialPoints( cbind( 
        sample(seq(bbox(layer)[1,1],bbox(layer)[1,2],length.out=100),nlocs),
        sample(seq(bbox(layer)[2,1],bbox(layer)[2,2],length.out=100),nlocs)
        ) )

nonmissing <- which(!is.na(values(layer)))

ndist <- 10

neighborhoods <- get.neighborhoods( ndist, locs, nonmissing, layer, numcores )
boundaries <- get.boundaries( neighborhoods, nonmissing, layer, numcores )
nonoverlapping <- which.nonoverlapping( neighborhoods )

nlayer <- nlayer.2 <- layer
values(nlayer)[nonmissing] <- ( ! seq_along(nonmissing) %in% unlist(neighborhoods) )
values(nlayer.2)[nonmissing] <- ( ! seq_along(nonmissing) %in% unlist(neighborhoods[nonoverlapping]) )
blayer <- blayer.2 <- layer
values(blayer)[nonmissing] <- ( ! seq_along(nonmissing) %in% unlist(boundaries) )
values(blayer.2)[nonmissing] <- ( ! seq_along(nonmissing) %in% unlist(boundaries[nonoverlapping]) )

layout(t(1:2))
plot(layer)
points(locs,pch="*",cex=2)
plot(nlayer,add=TRUE,col=adjustcolor(c("black",NA),.5))
plot(blayer,add=TRUE,col=adjustcolor(c("red",NA),.5))
plot(layer)
points(locs,pch="*",cex=2)
plot(nlayer.2,add=TRUE,col=adjustcolor(c("black",NA),.5))
plot(blayer.2,add=TRUE,col=adjustcolor(c("red",NA),.5))

layers <- sapply( layer.list, function (ll) { values(ll)[nonmissing] } )

ij <- adjacent(layer,cells=nonmissing,target=nonmissing,directions=4,pairs=TRUE,sorted=TRUE) # to and from cells both loc
ij <- ij[,2:1]
stopifnot( all(ij[,1] != ij[,2]) ) ## NO DIAGONAL

G <- sparseMatrix( i=match(ij[,1],nonmissing), j=match(ij[,2],nonmissing), x=1.0 )
Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )
stopifnot(nrow(layers)==nrow(G))

transfn <- exp
valfn <- function (gamma) { ( rowSums( layers * gamma[col(layers)], na.rm=TRUE ) ) }
ndelta <- ngamma <- nlayers
update.G <- function(params) {
    beta <- params[1]
    gamma <- params[1+(1:ngamma)]
    delta <- params[1+ngamma+(1:ndelta)]
    return( exp(beta) * transfn(valfn(gamma))[G@i+1L] * transfn( valfn(delta)[G@i+1L] + valfn(delta)[Gjj] ) )
    # return( exp(beta) * transfn(valfn(gamma))[G@i+1L] * transfn( abs(valfn(delta)[G@i+1L] - valfn(delta)[Gjj]) ) )
}

true.params <- c(3*runif(1),8*(runif(2*nlayers)-0.5))
names(true.params) <- c("beta",paste("gamma",1:nlayers,sep=''),paste("delta",1:nlayers,sep=''))

G@x <- update.G(true.params)

true.hts <- hitting.analytic( neighborhoods, G-diag(rowSums(G)), numcores=1 )

###
# check parameter inference
init.params <- true.params

# Massage the numerics.
zeros <- unlist(neighborhoods) + rep((seq_along(neighborhoods)-1)*nrow(true.hts),sapply(neighborhoods,length))
scaling <- 1 # sqrt(nrow(G) * length(locs))
sc.one <- 1/scaling
hts <- true.hts
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

update.aux(init.params,environment(),check=FALSE)
weightings <- ifelse( rowMeans(hts) < quantile(hts,.5), dG, 0 )  # indexes locations; note may overlap with zeros
nomitted <- sum( weightings[row(hts)[zeros]] )

L <- function (params) {
    update.aux(params,parent.env(environment()))
    ans <- ( sum( weightings*rowSums((GH+sc.one)^2) ) - (nomitted)*sc.one^2 )
    if (!is.finite(ans)) { browser() }
    return(ans)
}
dL <- function (params) {
    update.aux(params,parent.env(environment()))
    bgrad <- ( 2 )* sum( weightings * rowSums(GH * (GH+sc.one)) )
    ggrads <- sapply( 1:ncol(layers), function (kk) {
            2 * sum( weightings * rowSums( (layers[,kk] * GH) * (GH+sc.one)) )
        } )
    dgrads <- ggrads + sapply( 1:ncol(layers), function (kk) {
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

## NOTE that this won't look good at local minima (where the derivative is zero) i.e. at true.params
gcheck <- function (params=jitter(init.params),eps=1e-8,dp=eps*dirn,dirn=runif(length(params))) {
    # check gradient  (VERY STEEP ?!?!?!)
    L0 <- L(params)
    dL0 <- dL(params)
    L1 <- L(params+dp)
    dL1 <- dL(params+dp)
    c( L0, L1-L0, sum(dp*dL0), sum(dp*dL1) )
}

results <- list(par=true.params)

# check answer
layout(matrix(1:(2*(length(init.params)+1)),nrow=2,byrow=TRUE))
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


parscale <- c( abs(init.params[1]/10), rep(0.1,length(init.params)-1) )
# results <- optim( par=init.params, fn=L, gr=dL, control=list(parscale=parscale,fnscale=max(1,abs(L(init.params))/10)), method="BFGS" )


##
# The integral equation

IL <- function (params) {
    # integral.hts is the mean hitting time of each neighborhood to its boundary,
    #  plus the mean hitting time to each neighborhood (including itself)
    update.aux(params,parent.env(environment()))
    hitting.probs <- get.hitting.probs( G, dG, neighborhoods[nonoverlapping], boundaries[nonoverlapping], numcores=numcores )
    hitting.times <- get.hitting.times( G, dG, neighborhoods[nonoverlapping], boundaries[nonoverlapping], numcores=numcores )
    integral.hts <- do.call( rbind, mclapply( seq_along(neighborhoods[nonoverlapping]), function (k) {
            ihs <- hitting.times[[k]] + hitting.probs[[k]] %*% hts[boundaries[[nonoverlapping[k]]],nonoverlapping] 
            ihs[,k] <- 0
            return(ihs)
        }, mc.cores=numcores ) )
    ans <- sum( ( hts[unlist(neighborhoods[nonoverlapping]),nonoverlapping] / integral.hts - 1 )^2, na.rm=TRUE )
    # ans <- sum( ( hts[unlist(neighborhoods[nonoverlapping]),nonoverlapping] - integral.hts )^2, na.rm=TRUE )
    if (!is.finite(ans)) { browser() }
    return(ans)
}

# should be zero
IL(true.params)

# compare steepness
layout(matrix(1:((2*length(init.params)+1)),nrow=2))
fac <- 3
for (k in seq_along(init.params)) {
    parvals <- seq( init.params[k]-fac, init.params[k]+fac, length.out=20 )
    for (fn in list(L,IL)) {
        Lvals <- sapply(parvals, function (x) fn(ifelse(seq_along(init.params)==k,x,init.params)) )
        yrange <- range(1+Lvals,fn(init.params)+1)
        plot( parvals, 1+Lvals, ylim=yrange, main=names(init.params)[k], log='y' )
        abline(v=init.params[k])
        abline(h=fn(init.params))
    }
}


