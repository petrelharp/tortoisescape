#!/usr/bin/Rscript

source("resistance-fns.R")
require(raster)

layer.prefix <- c("geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_")
# layers <- c("annual_precip","barren_30","bd_ss2_st_30","eastness_30","lat_gcs_30","lon_gcs_30","slope30")
layer.names <- selfname( c("annual_precip","barren_30","lat_gcs_30") )
gamma <- c(.8, .001, .01)
delta <- c(.005,.4,.03)

# get info out of one layer
onelayer <- raster(paste(layer.prefix,layer.names[1],sep=''))
n <- dim(onelayer)[2]; m <- dim(onelayer)[1]
nmap <- matrix(1:n*m,nrow=n)
all.locs <- cbind( i=as.vector(row(nmap)), j=as.vector(col(nmap)) )

layers <- sapply(layer.names, function (ll) {
            rast <- raster(paste(layer.prefix,ll,sep=''))
            # note this is ROW-ORDERED
            # so to plot do:  dim(x) <- dim(rast)[2:1]; image(x)
            vrast <- scale( values(rast) )
            # FILL IN NAs AS ZEROS (for testing purposes)
            vrast[is.na(vrast)] <- 0
            return(vrast)
        } )
gc()
edges <- rBind( 
        grid.sum(n,m,direction=c(+1,0)) %*% layers, # rights
        grid.sum(n,m,direction=c(-1,0)) %*% layers, # lefts
        grid.sum(n,m,direction=c(0,+1)) %*% layers, # downs
        grid.sum(n,m,direction=c(0,-1)) %*% layers  # ups
    )

transfn <- exp
valfn <- function (gamma) { ( rowSums( sweep( layers, 2, gamma, "*" ), na.rm=TRUE ) ) }

update.G <- function(params) {
    gamma <- params[1:length(delta)]
    delta <- params[length(delta)+(1:length(gamma))]
    G <- grid.adjacency(n,m,diag=FALSE,symmetric=FALSE)
    stopifnot( length(G@x) == nrow(edges) )
    dp <- diff(G@p)
    jj <- rep(seq_along(dp),dp)
    G@x <- transfn(valfn(gamma))[G@i+1L] * transfn( valfn(delta)[G@i+1L] + valfn(delta)[jj] )
    diag(G) <- (-1) * rowSums(G)
    return(G)
}

G <- update.G(c(gamma,delta))
range(G@x)

# sampling locations
nsamps <- 10
locs.ij <- cbind( i=sample.int(n,nsamps)-1L, j=sample.int(m,nsamps)-1L )
locs <- ij.to.k(locs.ij,n,m)

# analytical mean hitting times
true.hts <- hitting.analytic(locs,G)  # warning, this could take a while (10s for n=100 and nsamps=20)

stopifnot(all(true.hts>=0))

# Observed hitting times
pairwise.hts <- true.hts[locs,]

if (FALSE) {
    # look at these
    for (kk in 1:ncol(true.hts)) {
        plot.ht( true.hts[,kk], dims=c(n,m) )
        if (is.null(locator(n=1))) { break; }
    }
}

###
# inference
g0 <- runif(length(gamma))/10
d0 <- runif(length(delta))/10

hts <- true.hts
zeros <- locs + (0:(length(locs)-1))*nrow(hts)

# Massage the numerics.
scaling <- sqrt(nrow(G) * nsamps)
sc.hts <- sweep(hts,2,colMeans(hts)*(nrow(hts)/(nrow(hts)-1)),"-")/scaling
sc.hts[zeros] <- 0
sc.one <- 1/scaling

#
GH <- G %*% sc.hts
GH[zeros] <- 0

L <- function (params) {
    if (any(params != get("params",parent.env(environment()) ) ) ) { 
        assign("params", params, parent.env(environment()) )
        assign( "G", update.G(params), parent.env(environment()) )
        GH <- G %*% sc.hts
        GH[zeros] <- 0
        assign("GH", GH, parent.env(environment()) )
    }
    return( sum( (GH+sc.one)^2 ) - length(zeros)*sc.one^2 )
}
dL <- function (params) {
    if (any(params != get("params", parent.env(environment()) ) ) ) { 
        assign("params", params, parent.env(environment()) )
        assign( "G", update.G(params), parent.env(environment()) )
        GH <- G %*% sc.hts
        GH[zeros] <- 0
        assign("GH", GH, parent.env(environment()) )
    }
    ggrads <- sapply( 1:length(gamma), function (kk) {
            2 * sum( layers[,kk] * GH * (GH+sc.one) )
        } )
    dgrads <- ggrads + sapply( 1:length(delta), function (kk) {
            GLH <- G %*% ( layers[,kk] * sc.hts )
            2 * sum( GLH * (GH+sc.one)  )
        } )
    return( c(ggrads, dgrads) )
}
environment(L) <- environment(dL) <- fun.env <- list2env( list(
                G=G,
                params=c(gamma,delta),
                GH=GH), 
        parent=environment() )

L(c(gamma,delta))
dL(c(gamma,delta))

L(c(gamma,delta)+.01)
dL(c(gamma,delta)+.01)

results <- optim( par=c(gamma,delta)+.1*runif(length(gamma)+length(delta),min=-1,max=1), fn=L, gr=dL, control=list(parscale=abs(c(gamma,delta))/200,trace=5), method="CG" )

G.est <- update.G(results$par)
est.hts <- hitting.analytic(locs,G.est)

range( (G.est%*%true.hts)[-zeros] )
range( est.hts - true.hts )


if (FALSE) {
    # look at these
    layout(1:2)
    for (kk in 1:ncol(true.hts)) {
        plot.ht( true.hts[,kk], dims=c(n,m) )
        plot.ht( est.hts[,kk], dims=c(n,m) )
        if (is.null(locator(n=1))) { break; }
    }
}

##
# check things
if (FALSE) {

    # check the gradient
    check.deriv <- function (p0,eps) {
        tmp <- c( L(p0+eps), L(p0-eps), L(p0) )
        c( analytic=sum(dL(p0)*eps), est1=(tmp[1]-tmp[3]), est2=(tmp[3]-tmp[2]) ) 
    }
    p0 <- c(gamma,delta)+.02
    for (k in seq_along(c(gamma,delta))) {
        eps <- 1e-6 * ifelse( k==seq_along(c(gamma,delta)), 1, 0 )
        print( check.deriv(p0,eps) )
    }
    eps <- runif( length(p0) ) * 1e-6
    check.deriv(p0,eps)

}
