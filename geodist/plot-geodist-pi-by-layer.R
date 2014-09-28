#!/usr/bin/Rscript

if (!interactive()) {
    subdir <- commandArgs(TRUE)[1]
    layer.name <- commandArgs(TRUE)[2]
} else {
    subdir <- "100x"
    layer.name <- "slope_30"
}


require(raster)
require(colorspace)
require(MASS)

# georaphic distance
torts <- read.csv("../1st_180_torts.csv",header=TRUE,stringsAsFactors=FALSE)
use.torts <- torts$EM_Tort_ID
nind <- nrow(torts)
tort.dist.table <- read.table("../1st180_pairwise_distances_sorted_redundancy_removed.txt",header=TRUE)
tort.dists <- numeric(nind^2); dim(tort.dists) <- c(nind,nind)
tort.dists[ cbind( match(tort.dist.table$etort1,torts$EM_Tort_ID), match(tort.dist.table$etort2,torts$EM_Tort_ID) ) ] <- tort.dist.table$DISTANCE
tort.dists <- tort.dists + t(tort.dists)

# environmental distance
edists <- read.table(paste(subdir,"/",layer.name,"-edist.tsv",sep=''),sep='\t',header=TRUE)

# coverages
coverage <- read.csv("../coverage_info.csv")
torts$coverage <- coverage$sequence_yield_megabases[match(torts$EM_Tort_ID,coverage$individual)]

# genetic distance
pimat.vals <- scan("../pairwisePi/alleleCounts_1millionloci.pwp") # has UPPER with diagonal
pimat <- numeric(nind^2)
dim(pimat) <- c(nind,nind)
pimat[upper.tri(pimat,diag=TRUE)] <- pimat.vals
pimat[lower.tri(pimat,diag=FALSE)] <- t(pimat)[lower.tri(pimat,diag=FALSE)]
dimnames(pimat) <- list( torts$Tort_EM_Id, torts$Tort_EM_Id )

# combine them
dists <- merge( tort.dist.table, edists, by.x=c("etort1","etort2"), by.y=c("tort1","tort2") )
dists$pi <- pimat[ cbind(match(dists$etort1,torts$EM_Tort_ID),match(dists$etort2,torts$EM_Tort_ID)) ]

# remove missing one
dists <- droplevels(na.omit(dists))
use.torts <- levels(torts$EM_Tort_ID)[levels(torts$EM_Tort_ID) %in% c(levels(dists$etort1),levels(dists$etort2))]
dists$etort1 <- factor( as.character(dists$etort1), use.torts )
dists$etort2 <- factor( as.character(dists$etort2), use.torts )

ntorts <- length(use.torts)
rm(nind)

##
# scale
dists$DISTANCE <- dists$DISTANCE/1000

pr <- function (params) {
    A <- params[1]; B <- params[2]; Const <- params[3]; C <- params[3+(1:ntorts)]
    Z <- with(dists, C[etort1] * C[etort2] * ( A * DISTANCE + B * edist + Const ) )
    plot( pi ~ Z, data=dists )
    abline(0,1,col='red')
    with( dists, mean( ( C[etort1] * C[etort2] * ( A * DISTANCE + B * edist + Const ) - pi )^2 ) )
}

base.lm <- lm( pi ~ DISTANCE + edist, data=dists )
init.params <- c( A=coef(base.lm)[2], B=coef(base.lm)[3], Const=coef(base.lm)[1], C=rep(1,length(use.torts)) )

pr(init.params)

## want to do this: but nls won't
# fit <- nls( pi ~  C[etort1] * C[etort2] * ( A * DISTANCE + B * edist + Const ), data=dists,  ...
#
# NOTE: constrain C to have C[1] = 1 or else nonidentifiable


params <- init.params

for (k in 1:160) {
    cat(k,"\n")
    old.params <- params

    # Given C, fit A, B, Const:
    C <- c(1,params[3+(1:ntorts)])
    AB.lm <- lm( I(pi/(C[etort1]*C[etort2])) ~ DISTANCE + edist, data=dists )
    params[1:3] <- coef(AB.lm)[c(2,3,1)]

    pr(params)

    # given A, B, Const, fit C
    C.model.matrix1 <- model.matrix( pi ~ etort1, data=dists )
    colnames(C.model.matrix1) <- gsub("^etort[12]","",colnames(C.model.matrix1))
    Cmm1 <- C.model.matrix1[,match(levels(dists$etort1),colnames(C.model.matrix1))]
    Cmm1[is.na(Cmm1)] <- 0
    C.model.matrix2 <- model.matrix( pi ~ etort2, data=dists )
    colnames(C.model.matrix2) <- gsub("^etort[12]","",colnames(C.model.matrix2))
    Cmm2 <- C.model.matrix1[,match(levels(dists$etort1),colnames(C.model.matrix2))]
    Cmm2[is.na(Cmm2)] <- 0
    stopifnot(all(levels(dists$etort1)==levels(dists$etort2)))
    C.model.matrix <- Cmm1 + Cmm2

    # minimize ( z - A x )^2 by solving  A^T  z = A^T A x
    A <- params[1]; B <- params[2]; Const <- params[3]
    lr.pi <- with(dists, log(pi / (A*DISTANCE+B*edist+Const)) )

    log.C.lm <- ginv(crossprod(C.model.matrix)) %*% crossprod(C.model.matrix,lr.pi)
    C.lm <- exp(log.C.lm)
    params[-(1:3)] <- C.lm

    pr(params)

    if (max((old.params-params)*c(1e6,1e6,1,rep(1,length(C)-1)))<1e-8) { break }
}

write( c("DISTANCE","edist",use.torts), file=paste(subdir,"/",layer.name,"-inferred-geodist-params.tsv",sep='' ), ncolumns=length(params) )
write( params, file=paste(subdir,"/",layer.name,"-inferred-geodist-params.tsv",sep='' ), ncolumns=length(params), append=TRUE )

# look at results
A <- params[1]; B <- params[2]; Const <- params[3]
C <- c(1,params[3+(1:ntorts)])
dists$npi <- dists$pi / ( C[dists$etort1] * C[dists$etort2] )

layout(t(1:2))
with(dists, plot( DISTANCE, pi ) )
with(dists, plot( A*DISTANCE + B*edist, npi ) )

layout(t(1:2))
with(dists, {
        plot( A*DISTANCE + B*edist, npi, type='n' )
        text( A*DISTANCE + B*edist, npi, labels=gsub("etort.","",dists$etort1) )
    })
with(dists, {
        plot( A*DISTANCE + B*edist, npi, type='n' )
        text( A*DISTANCE + B*edist, npi, labels=gsub("etort.","",dists$etort2) )
    })

##
# look at wierdos
get.neighbors <- function (inds) {
    sapply( inds, function (w) {
        sdists <- subset(dists,etort1==w | etort2==w)
        kmin <- which.min( sdists$DISTANCE )
        with( sdists, if (etort1[kmin]==w) { etort2[kmin] } else { etort1[kmin] } )
    } )
}

wierd.counts <- with( subset(dists, npi < .33), table(etort1)+table(etort2) )
wierdos <- names(wierd.counts)[wierd.counts>10]
wierd.neighbors <- get.neighbors(wierdos)
normals <- sample(use.torts,length(wierdos))
normal.neighbors <- get.neighbors(normals)

##
# plotting whatnot
layer <- raster("../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_dem_30")
load("../tort.coords.rasterGCS.Robj")
coord.names <- rownames(tort.coords.rasterGCS@coords)

pair.segments <- function (ind1,ind2) {
    wd <- subset(dists, etort1==ind1 | etort2==ind1)
    wd$etort <- with(wd, levels(etort1)[ifelse( etort1==ind1, etort2, etort1 )] )
    nd <- subset(dists, etort1==ind2 | etort2==ind2)
    nd$etort <- with(nd, levels(etort1)[ifelse( etort1==ind2, etort2, etort1 )] )
    wd <- wd[!wd$etort==ind2,]
    nd <- nd[!nd$etort==ind1,]
    nd <- nd[ match(wd$etort,nd$etort), ]
    stopifnot(all( wd$etort==nd$etort ) )
    with( wd, points( DISTANCE, pi, pch=20, col='red' ) )
    with( nd, points( DISTANCE, pi, pch=20, col='green' ) )
    segments( x0=wd$DISTANCE, x1=nd$DISTANCE, y0=wd$pi, y1=nd$pi, col='purple' )
}

pdf(file=paste(subdir,"/",layer.name,"-wierdos.pdf",sep=''), width=8, height=8, pointsize=10 )
layout(matrix(1:4,nrow=2))
for (k in seq_along(wierdos) ) {
    plot(layer,main="wierd")
    text(tort.coords.rasterGCS, gsub("etort.","",coord.names), 
        col=ifelse( coord.names%in%wierdos[k], "red", ifelse( coord.names%in%wierd.neighbors[k], "green", "black") ) )
    with(dists, plot( DISTANCE, pi, pch=20, cex=0.5, col=adjustcolor("black",0.5) ) )
    pair.segments(wierdos[k],wierd.neighbors[k])
    #
    plot(layer,main="normal")
    text(tort.coords.rasterGCS, gsub("etort.","",coord.names), 
        col=ifelse( coord.names%in%normals[k], "red", ifelse( coord.names%in%normal.neighbors[k], "green", "black") ) )
    with(dists, plot( DISTANCE, pi, pch=20, cex=0.5, col=adjustcolor("black",0.5) ) )
    pair.segments(normals[k],normal.neighbors[k])
}
dev.off()


if (FALSE) {

## numerical issues!!

    dists$norm.pi <- dists$pi - mean(dists$pi)

H <- function (params) {
    A <- params[1]
    B <- params[2]
    Const <- params[3]
    C <- c(1, params[3+(1:(ntorts-1))] )
    with( dists, sum( ( norm.pi - C[etort1] * C[etort2] * ( A * DISTANCE + B * edist + Const ) )^2 ) )
}
dH <- function (params) {
    A <- params[1]
    B <- params[2]
    Const <- params[3]
    C <- c(1, params[3+(1:(ntorts-1))])
    CC <- with(dists, C[etort1] * C[etort2] )
    Z <- with(dists, ( norm.pi - CC * ( A * DISTANCE + B * edist + Const ) ) )
    dC <- with(dists, cbind( 
            tapply( C[etort2] * ( A * DISTANCE + B * edist + Const ) * Z, etort1, sum, na.rm=TRUE ),
            tapply( C[etort1] * ( A * DISTANCE + B * edist + Const ) * Z, etort2, sum, na.rm=TRUE ) 
        )[-1,] )
    with( dists, -2 * c( 
            sum( DISTANCE * CC * Z ), # A
            sum( edist * CC * Z ), # B
            sum( CC * Z ), # Const
            rowSums(dC,na.rm=TRUE) #C
            )
        )
}

base.lm <- lm( norm.pi ~ DISTANCE + edist, data=dists )
init.params <- c( A=coef(base.lm)[2], B=coef(base.lm)[3], Const=coef(base.lm)[1], C=rep(1,ntorts-1) ) 
parscale <- c( 1e3, 3e2, 1, rep(.01,ntorts-1) )

grad.check <- t( sapply( seq_along(init.params), function (k) {
        eps <- rep(0,length(init.params))
        eps[k] <- 1e-4 
        ( c(
            H(init.params),
            H(init.params+eps) - H(init.params),
            sum( dH(init.params) * eps )
        ) )
    } ) )


soln <- optim( par=init.params, fn=H, gr=dH, method="BFGS", control=list(parscale=parscale) )


}
