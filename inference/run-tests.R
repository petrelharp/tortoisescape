#!/usr/bin/Rscript

require(raster)
source("resistance-fns.R")

biglayer <- raster(nrows=2^10,ncols=2^10,xmn=0,ymn=0,xmx=100,ymx=100)
values(biglayer) <- rnorm(length(biglayer))
nlayers <- 2
layer.list <- lapply( 1:nlayers, function (k) {
        values(biglayer) <- rnorm(length(biglayer))
        aggregate( biglayer, fact=2^7, fun=mean )
    } )
layer <- layer.list[[1]]

nlocs <- 4
locs <- SpatialPoints( cbind( 
        sample(seq(bbox(layer)[1,1],bbox(layer)[1,2],length.out=100),nlocs),
        sample(seq(bbox(layer)[2,1],bbox(layer)[2,2],length.out=100),nlocs)
        ) )

nonmissing <- seq_along(values(layer))

ndist <- 3

neighborhoods <- lapply( seq_along(locs) , function (k) {
        d_tort <- distanceFromPoints( layer, locs[k] )
        match( Which( d_tort <= max(ndist,minValue(d_tort)), cells=TRUE, na.rm=TRUE ), nonmissing )
    } )
nlayer <- layer
values(nlayer) <- ( ! seq_along(values(layer)) %in% unlist(neighborhoods) )

plot(layer)
plot(nlayer,add=TRUE,col=adjustcolor(c("black",NA),.5))
points(locs,pch="*",cex=2)

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
    return( beta * transfn(valfn(gamma))[G@i+1L] * transfn( valfn(delta)[G@i+1L] + valfn(delta)[Gjj] ) )
}

true.params <- 2*(runif(1+2*nlayers)-0.5)
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
    bgrad <- ( 2 / params[1] )* sum( weightings * rowSums(GH * (GH+sc.one)) )
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

gcheck <- function (params=init.params,eps=1e-8,dp=eps*dirn,dirn=runif(length(params))) {
    # check gradient  (VERY STEEP ?!?!?!)
    L0 <- L(params)
    dL0 <- dL(params)
    L1 <- L(params+dp)
    dL1 <- dL(params+dp)
    c( L0, L1-L0, sum(dp*dL0), sum(dp*dL1) )
}

results <- list(par=true.params)

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



parscale <- c( abs(init.params[1]/10), rep(0.1,length(init.params)-1) )
results <- optim( par=init.params, fn=L, gr=dL, control=list(parscale=parscale,fnscale=max(1,abs(L(init.params))/10)), method="BFGS" )
