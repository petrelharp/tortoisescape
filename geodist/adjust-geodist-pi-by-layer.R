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

# geographic distance
torts <- read.csv("../tort_180_info/1st_180_torts.csv",header=TRUE,stringsAsFactors=FALSE)
nind <- nrow(torts)
tort.dist.table <- read.table("../tort_180_info/1st180_pairwise_distances_sorted_redundancy_removed.txt",header=TRUE)
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
use.torts <- torts$EM_Tort_ID[torts$EM_Tort_ID %in% c(levels(dists$etort1),levels(dists$etort2))]
dists$etort1 <- factor( as.character(dists$etort1), use.torts )
dists$etort2 <- factor( as.character(dists$etort2), use.torts )

ntorts <- length(use.torts)
rm(nind)

##
# scale
dists$DISTANCE <- dists$DISTANCE/1000


## want to do this: but nls won't
# fit <- nls( pi ~  C[etort1] * C[etort2] * ( A * DISTANCE + B * edist + Const ), data=dists,  ...
#
##
# resolve nonidentifiability by setting Const = mean(pi) .
mean.pi <- mean(dists$pi)

pr <- function (params) {
    A <- params[1]; B <- params[2]; C <- params[2+(1:ntorts)]
    Z <- with(dists, C[etort1] * C[etort2] * ( A * DISTANCE + B * edist + mean.pi ) )
    with( dists, plot( Z,  pi ) )
    abline(0,1,col='red')
    cat(with( dists, mean( ( C[etort1] * C[etort2] * ( A * DISTANCE + B * edist + mean.pi ) - pi )^2 ) ),"\n")
}



# setup for given A, fit C
C.model.matrix <- sapply( 1:ntorts, function (k) {
        ( dists$etort1 == use.torts[k] ) + ( dists$etort2 == use.torts[k] )
    } )
colnames(C.model.matrix) <- use.torts

base.lm <- lm( pi-mean.pi ~ DISTANCE+edist+0, data=dists )
init.params <- c( A=coef(base.lm)[1], B=coef(base.lm)[2], C=rep(1,length(use.torts)) )
names(init.params) <- c("DISTANCE","edist",use.torts)

pr(init.params)

params <- init.params

max.k <- 100
for (k in 1:max.k) {
    cat(k,"\n")
    old.params <- params

    # Given C, fit A, B
    C <- params[2+(1:ntorts)]
    AB.lm <- lm( I(pi/(C[etort1]*C[etort2])-mean.pi) ~ DISTANCE + edist + 0, data=dists )
    params[1:2] <- coef(AB.lm)[1:2]

    pr(params)

    # minimize ( z - A x )^2 by solving  A^T  z = A^T A x
    A <- params[1]; B <- params[2] 
    lr.pi <- with(dists, log(pi / (A*DISTANCE+B*edist+mean.pi)) )

    log.C.lm <- ginv(crossprod(C.model.matrix)) %*% crossprod(C.model.matrix,lr.pi)
    params[-(1:2)] <- exp(log.C.lm)

    pr(params)

    if (max(abs(old.params-params)*c(1e6,1,rep(1,length(C))))<1e-4) { break }
}
if (k==max.k) { warning("Didn't converge.") }


outfile <- paste(subdir,"/",layer.name,"-inferred-geodist-params.tsv",sep='' )
write( c("DISTANCE","edist",use.torts), file=outfile, ncolumns=length(params) )
write( params, file=outfile, ncolumns=length(params), append=TRUE )

# look at results
A <- params[1]; B <- params[2]
C <- params[2+(1:ntorts)]
dists$npi <- dists$pi / ( C[dists$etort1] * C[dists$etort2] )

write.csv( dists, file=paste(subdir,"/",layer.name,"-pairwise-normalized-pi.csv", sep='' ), row.names=FALSE )

pdf( file=paste(subdir,"/",layer.name,"-adjust-pi-results.pdf", sep='' ), width=8,height=4,pointsize=10)

layout(t(1:2))
with(dists, plot( DISTANCE, pi, pch=20, cex=0.5 ) )
with(dists, plot( DISTANCE, npi, pch=20, cex=0.5 ) )

##
# plotting whatnot
layer <- raster("../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_dem_30")
load("../tort_180_info/tort.coords.rasterGCS.Robj")
coord.names <- rownames(tort.coords.rasterGCS@coords)

plot(layer,main="adjustment factor")
text(tort.coords.rasterGCS, gsub("etort.","",coord.names),cex=.5)
plot(layer,main="adjustment factor")
points(tort.coords.rasterGCS, cex=(C[match(coord.names,names(C))]-.75)/.10 )

dev.off()

