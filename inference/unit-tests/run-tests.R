#!/usr/bin/Rscript
require(numDeriv)

# will compare results to those previously computed with the same seed
set.seed(12345)
check.file <- "run-tests-saved-1234.RData"
check.objects <- c("biglayer","layer.list","SP.locs","locs","nonmissing","neighborhoods","boundaries","nonoverlapping","G","Gjj","true.hts","sub.true.hts","gcheck.1","gcheck.2","hcheck.1","hcheck.2","hts.0","hts.1","hts.2","hts.3","hts.4","hts.5", "hs", "hs.checks", "hc")
# do this at the end
check.it <- function () { 
    if (do.check) {
        checks <- sapply( selfname(check.objects), function (check.me) {
                all.equal( get(check.me), get(check.me,envir=saved.env) ) } )
        log.checks <- ( as.character(checks) == "TRUE" )
        if ( any( ! log.checks ) ) {
            stop( paste( paste( check.objects[!log.checks], checks[!log.checks], sep=" : " ), collapse="\n" ) )
        } else { cat("\nrun-tests.R:\n  Everything checks out.\n") }
    } else  {
        save( list=check.objects, file=check.file )
        cat("Saved output in ", check.file, "\n")
    }
}

do.check <- file.exists(check.file)
if (do.check) {
    saved.env <- new.env(parent=emptyenv())
    saved.objects <- load( file=check.file, env=saved.env )
    if( !all( check.objects %in% saved.objects ) ) { stop(paste("Missing", paste(setdiff(check.objects,saved.objects),collapse=", "), "from", check.file) ) }
}

require(raster)
source("../resistance-fns.R")
require(parallel); numcores <- getcores()

# set up some layers
biglayer <- raster(nrows=2^10,ncols=2^10,xmn=0,ymn=0,xmx=100,ymx=100)
values(biglayer) <- rnorm(length(biglayer))
nlayers <- 2
fact <- 2^5
layer.list <- lapply( 1:nlayers, function (k) {
        values(biglayer) <- rnorm(length(biglayer))
        aggregate( biglayer, fact=2^5, fun=mean, na.rm=FALSE )
    } )
# include some missing values
navals <- sample(length(values(layer.list[[1]])),10)
for (k in 1:nlayers) { values(layer.list[[k]])[navals] <- NA }
layer <- layer.list[[1]]

# sampling locations
nlocs <- 40
SP.locs <- SpatialPoints( cbind( 
        sample(seq(bbox(layer)[1,1],bbox(layer)[1,2],length.out=100),nlocs),
        sample(seq(bbox(layer)[2,1],bbox(layer)[2,2],length.out=100),nlocs)
        ) )

nonmissing <- which(!is.na(values(layer)))

orig.locs <- cellFromXY( layer, SP.locs )
locs <- match(orig.locs,nonmissing)

ndist <- 10

neighborhoods <- get.neighborhoods( ndist, SP.locs, nonmissing, layer, numcores )
boundaries <- get.boundaries( neighborhoods, nonmissing, layer, numcores )
nonoverlapping <- which.nonoverlapping( neighborhoods )

# neighborhoods make sense?
stopifnot( all( sapply( seq_along(locs), function (k) { locs[k] %in% neighborhoods[[k]] } ) ) )

nlayer <- nlayer.2 <- layer
values(nlayer)[nonmissing] <- ( ! seq_along(nonmissing) %in% unlist(neighborhoods) )
values(nlayer.2)[nonmissing] <- ( ! seq_along(nonmissing) %in% unlist(neighborhoods[nonoverlapping]) )
blayer <- blayer.2 <- layer
values(blayer)[nonmissing] <- ( ! seq_along(nonmissing) %in% unlist(boundaries) )
values(blayer.2)[nonmissing] <- ( ! seq_along(nonmissing) %in% unlist(boundaries[nonoverlapping]) )

if (interactive()) {  # visual check of neighborhoods
    layout(t(1:2))
    plot(layer)
    points(SP.locs,pch="*",cex=2)
    plot(nlayer,add=TRUE,col=adjustcolor(c("black",NA),.5))
    plot(blayer,add=TRUE,col=adjustcolor(c("red",NA),.5))
    plot(layer)
    points(SP.locs,pch="*",cex=2)
    plot(nlayer.2,add=TRUE,col=adjustcolor(c("black",NA),.5))
    plot(blayer.2,add=TRUE,col=adjustcolor(c("red",NA),.5))
}

layers <- sapply( layer.list, function (ll) { values(ll)[nonmissing] } )

ij <- adjacent(layer,cells=nonmissing,target=nonmissing,directions=4,pairs=TRUE,sorted=TRUE) # to and from cells both loc
ij <- ij[,2:1]
stopifnot( all(ij[,1] != ij[,2]) ) ## NO DIAGONAL

G <- sparseMatrix( i=match(ij[,1],nonmissing), j=match(ij[,2],nonmissing), x=1.0 )
Gjj <- rep( seq.int(length(G@p)-1), diff(G@p) )
stopifnot(nrow(layers)==nrow(G))

transfn <- function (x) { 1/(1+exp(-x)) }
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

# analytically computed hitting times
true.hts <- hitting.analytic( neighborhoods, G-diag(rowSums(G)), numcores=1 )

# can compute subsets?
check.inds <- 2:4
sub.true.hts <- hitting.analytic( neighborhoods[check.inds], G-diag(rowSums(G)) )
stopifnot(all( abs( true.hts[,check.inds]-sub.true.hts ) < 1e-16 ) )

###
# check parameter inference
init.params <- true.params

# Massage the numerics.
zeros <- unlist(neighborhoods) + rep((seq_along(neighborhoods)-1)*nrow(true.hts),sapply(neighborhoods,length))
scaling <- 1 # sqrt(nrow(G) * length(locs))
sc.one <- 1/scaling
hts <- true.hts
hts[zeros] <- 0

# obejctive funciton and gradient for parameter inference
LdL <- params.logistic.setup(init.params,G,update.G,hts,zeros,sc.one,layers,transfn,valfn,ndelta,ngamma)
L <- LdL$L
dL <- LdL$dL

L(init.params)
dL(init.params)

# check that dL is the gradient of L
gcheck.1 <- gcheck(f=L, df=dL, params=init.params, eps=1e-8)
stopifnot( all( abs( gcheck.1[,3] - gcheck.1[,4] ) < 1e-8 ) )
gcheck.2 <- gcheck(f=L, df=dL, params=init.params+0.1, eps=1e-8)
stopifnot( all( abs( gcheck.2[,3] - gcheck.2[,4] ) < 1e-8 ) )

##
# check the hessian computation
LddL <- logistic.trust.setup(G,update.G,hts,zeros,sc.one,layers,transfn,valfn,ndelta,ngamma)
hcheck.1 <- hcheck(f=LddL, params=init.params, eps=1e-6)
hcheck.2 <- hcheck(f=LddL, params=init.params+0.1, eps=1e-6)

check.hcheck <- function (hc) {
    stopifnot( all( abs( hc["diff",] - hc["df0",] ) < 1e-6 ) && all( abs( hc["diff",] - hc["df1",] ) < 1e-6 ) )  # first order diff
    stopifnot( all( abs( hc["diff",] - hc["df2.0",] ) < 1e-9 ) && all( abs( hc["diff",] - hc["df2.1",] ) < 1e-9 ) )  # second order
    for (k in 1:(2*nlayers+1)) {
        stopifnot( all( abs( hc[6+(2*nlayers+1)+k,] - hc[6+k+2*(2*nlayers+1),] ) < 1e-6 ) && all( abs( hc[6+k+(2*nlayers+1),] - hc[6+k+3*(2*nlayers+1),] ) < 1e-6 ) ) # diff for grad
    }
}
check.hcheck(hcheck.1)
check.hcheck(hcheck.2)


###########
# check hitting time interpolation

obs.ht <- true.hts[locs,]

# these should be the same as they start from the truth
hts.0 <- interp.hitting( neighborhoods, G-diag(rowSums(G)), obs.ht, obs.locs=locs, alpha=0.0, numcores=numcores )
hts.1 <- interp.hitting( neighborhoods, G-diag(rowSums(G)), obs.ht, obs.locs=locs, alpha=1.0, numcores=numcores )

stopifnot(all(abs(hts.0-true.hts)<1e-7))
stopifnot(all(abs(hts.0-hts.1)<1e-7))

# and add noise
eps <- .01
noisy.ht <- true.hts[locs,] * exp( rnorm(length(obs.ht)) * eps )
hts.2 <- interp.hitting( neighborhoods, G-diag(rowSums(G)), noisy.ht, obs.locs=locs, alpha=0.0, numcores=numcores )
hts.3 <- interp.hitting( neighborhoods, G-diag(rowSums(G)), noisy.ht, obs.locs=locs, alpha=1.0, numcores=numcores )

stopifnot( all.equal( hts.0, hts.2 ) )
stopifnot( all( abs( ((hts.3-true.hts)/true.hts)[true.hts>0] ) < 5*eps ) )

# and more noise
eps <- .1
noisy.ht <- true.hts[locs,] * exp( rnorm(length(obs.ht)) * eps )
hts.4 <- interp.hitting( neighborhoods, G-diag(rowSums(G)), noisy.ht, obs.locs=locs, alpha=0.0, numcores=numcores )
hts.5 <- interp.hitting( neighborhoods, G-diag(rowSums(G)), noisy.ht, obs.locs=locs, alpha=1.0, numcores=numcores )

stopifnot( all.equal( hts.0, hts.4 ) )
stopifnot( all( abs( ((hts.5-true.hts)/true.hts)[true.hts>0] ) < 5*eps ) )


###
# check derivative of hitting times wrt parameters


hs <- hitting.sensitivity(true.params, neighborhoods, G, update.G, layers, transfn, valfn, ndelta, ngamma, do.hessian=TRUE, numcores=numcores)
epsval <- 1e-6
hs.grad.checks <- lapply( seq_along(true.params), function (k) {
        eps <- epsval * ifelse( seq_along(true.params)==k, 1, 0 )
        G@x <- update.G(true.params+eps)
        d.true.hts <- hitting.analytic( neighborhoods, G-diag(rowSums(G)), numcores=numcores )
        return( ( d.true.hts - true.hts ) - epsval * hs$gradient[,,k] )
    } )
stopifnot(all(abs(unlist(hs.grad.checks)<1e-9)))

num.jacob <- jacobian( func=function(x) { G@x <- update.G(x); as.vector(hitting.analytic( neighborhoods, G-diag(rowSums(G)), numcores=numcores )) }, x=true.params )
stopifnot( all( abs( num.jacob - matrix(hs$gradient,ncol=length(true.params)) ) < 1e-5 ) )

# check hessian of hitting times wrt parameters
check.loc <- c( neighborhoods[[1]][1] + 20, 1 )
f <- function (params) {
    G@x <- update.G(params)
    hitting.analytic( neighborhoods, G-diag(rowSums(G)), numcores=numcores )[check.loc[1],check.loc[2]]
}

num.grad <- grad( func=f, x=true.params )
stopifnot( all( abs( num.grad - hs$gradient[check.loc[1],check.loc[2],] ) < 1e-5 ) ) 

num.hess <- hessian( func=f, x=true.params )
stopifnot( all( abs( num.hess - hs$hessian[check.loc[1],check.loc[2],,] ) < 1e-6 ) )


# repeat a layer to create collinearity
coll.layers <- cbind( layers, layers[,1] )
coll.params <- c( coll.params[1], coll.params[1+(1:nlayers)], 0, coll.params[1+nlayers+(1:nlayers)], 0 )
ngamma <- nlayers+1
ndelta <- nlayers+1
# this should now be degenerate
hc <- hitting.colinearity(params=coll.params, locs=neighborhoods, obs.locs=locs, G=G, update.G=update.G, layers=coll.layers, transfn=transfn, valfn=valfn, ndelta=ndelta, ngamma=ngamma, numcores=numcores)

cvec.1 <- numeric(length(coll.params))
cvec.1[2] <- 1
cvec.1[2+nlayers] <- (-1)
stopifnot( all.equal( as.numeric(cvec.1 %*% hc %*% cvec.1), 0 ) )

cvec.2 <- numeric(length(coll.params))
cvec.2[3+nlayers] <- 1
cvec.2[3+2*nlayers] <- (-1)
stopifnot( all.equal( as.numeric(cvec.2 %*% hc %*% cvec.2), 0 ) )


###
# check everything agrees with previously saved versions
check.it()
