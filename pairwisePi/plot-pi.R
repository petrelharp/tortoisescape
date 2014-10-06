#!/usr/bin/Rscript
require(rgdal)
require(maps)
require(maptools)
require(raster)

angsd.pimat.vals <- scan("350000.pwp")  # has LOWER triangle of entries, without diagonal
robust.pimat.vals <- scan("alleleCounts_1millionloci.pwp") # robust version, has UPPER with diagonal
unbiased1.pimat.vals <- scan("unbiased/1millionALLsites.pwp")  # as robust version
unbiased10.pimat.vals <- scan("unbiased/10millionALLsites.pwp")  # as robust version
nind <- 180

robust.pimat <- angsd.pimat <- unbiased1.pimat <- unbiased10.pimat <- numeric(nind^2)
dim(robust.pimat) <- dim(angsd.pimat) <- dim(unbiased1.pimat) <- dim(unbiased10.pimat) <- c(nind,nind)
angsd.pimat[lower.tri(angsd.pimat,diag=FALSE)] <- angsd.pimat.vals
angsd.pimat <- angsd.pimat + t(angsd.pimat)
robust.pimat[upper.tri(robust.pimat,diag=TRUE)] <- robust.pimat.vals
robust.pimat[lower.tri(robust.pimat,diag=FALSE)] <- t(robust.pimat)[lower.tri(robust.pimat,diag=FALSE)]
unbiased1.pimat[upper.tri(unbiased1.pimat,diag=TRUE)] <- unbiased1.pimat.vals
unbiased1.pimat[lower.tri(unbiased1.pimat,diag=FALSE)] <- t(unbiased1.pimat)[lower.tri(unbiased1.pimat,diag=FALSE)]
unbiased10.pimat[upper.tri(unbiased10.pimat,diag=TRUE)] <- unbiased10.pimat.vals
unbiased10.pimat[lower.tri(unbiased10.pimat,diag=FALSE)] <- t(unbiased10.pimat)[lower.tri(unbiased10.pimat,diag=FALSE)]

mean(angsd.pimat)
mean(robust.pimat)
mean(unbiased1.pimat)
mean(unbiased10.pimat)

torts <- read.csv("../1st_180_torts.csv",header=TRUE)
torts$EM_Tort_ID <- levels(torts$EM_Tort_ID)[as.numeric(torts$EM_Tort_ID)]
torts$EM_Tort_ID <- factor( torts$EM_Tort_ID , levels=torts$EM_Tort_ID  )

# pairwise distances
dists <- read.table("../1st180_pairwise_distances_sorted_redundancy_removed.txt",header=TRUE,stringsAsFactors=FALSE)
dists$etort1 <- factor( dists$etort1, levels=levels(torts$EM_Tort_ID) )
dists$etort2 <- factor( dists$etort2, levels=levels(torts$EM_Tort_ID) )
tort.dists <- numeric(nind^2); dim(tort.dists) <- c(nind,nind)
tort.dists[ cbind( match(dists$etort1,torts$EM_Tort_ID), match(dists$etort2,torts$EM_Tort_ID) ) ] <- dists$DISTANCE
tort.dists <- tort.dists + t(tort.dists)

# coverage
coverages <- read.csv("../coverage_info.csv")
torts$coverage <- coverages$sequence_yield_megabases

# elevation raster
elev.file <- "../geolayers/TIFF/10x/crop_resampled_masked_aggregated_10x_dem_30.gri"
elev <- raster(elev.file)
# and tortoise locs on it
load("../tort.coords.rasterGCS.Robj")
load("../county_lines.Robj")  # @gbradburd: how was this produced?

# make a table
dists$angsd <- angsd.pimat[ cbind( dists$etort1, dists$etort2 ) ]
dists$robust <- robust.pimat[ cbind( dists$etort1, dists$etort2 ) ]
dists$unbiased1 <- unbiased1.pimat[ cbind( dists$etort1, dists$etort2 ) ]
dists$unbiased10 <- unbiased10.pimat[ cbind( dists$etort1, dists$etort2 ) ]




# compare the methods
ut <- upper.tri(tort.dists,diag=FALSE) 
pdf(file="pi-methods-comparison.pdf",width=12,height=12,pointsize=10)
lcols <- adjustcolor(rainbow(nlevels(torts$Location_ID)),0.75)
pairs( dists[,-(1:2)], pch=20, cex=0.5, 
        upper.panel=function (x,y,...) points(x,y,col=lcols[torts$Location_ID][row(pimat)[ut]],...),
        lower.panel=function (x,y,...) points(x,y,col=lcols[torts$Location_ID][row(pimat)[ut]],...)
        )
dev.off()

stop('here')

###
# make plots for each method
for (meth in c("angsd","robust","unbiased1","unbiased10")) {

    pimat <- get(paste(meth,"pimat",sep='.'))

    ut <- upper.tri(pimat,diag=FALSE) 
    require(colorspace)
    tcols <- rainbow(240)
    lcols <- adjustcolor(rainbow(nlevels(torts$Location_ID)),0.75)

    ibd.list <- lapply( 1:nind, function (k) lm( pimat[k,] ~ tort.dists[k,] ) )
    ibd.coef <- sapply(ibd.list,coef)

    var_col <- function (varname) { 1 + outer( torts[[varname]], torts[[varname]], "!=" ) }

    torts$Adapter11 <- ( torts$Adapter == "BFIDT-011" )

    if (FALSE) {
        for (varname in c("Sex","Cleaned","Tissue.Box","Adapter11","Lib.prepped.","Sent.to.Berk.")) {
            plot( tort.dists[ut], pimat[ut], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5, col=var_col(varname)[ut], main=varname )
            if (is.null(locator(1))) { break }
        }
    }

    et110 <- which( torts$EM_Tort_ID == "etort-110" )
    et50 <- which( torts$EM_Tort_ID == "etort-50" )
    other11 <- which( torts$EM_Tort_ID != "etort-110" & ! torts$EM_Tort_ID == "etort-50" & torts$Adapter11 ) 

    pdf(file=paste(meth,"pairwise-pi.pdf",sep='-'), width=10, height=6, pointsize=10)

        plot( tort.dists[ut], pimat[ut], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5 )
        abline( coef( lm( pimat[ut] ~ tort.dists[ut] ) ) ) 

        plot( tort.dists[ut], pimat[ut], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5, col=lcols[torts$Location_ID][row(pimat)[ut]] ) # tcols[row(pimat)[ut]] )
        abline( coef( lm( pimat[ut] ~ tort.dists[ut] ) ) ) 
        legend("bottomright",pch=20,col=lcols,legend=levels(torts$Location_ID),cex=0.5)

        plot( tort.dists[ut], pimat[ut], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5, col=lcols[torts$Location_ID][col(pimat)[ut]] ) # tcols[col(pimat)[ut]] )
        abline( coef( lm( pimat[ut] ~ tort.dists[ut] ) ) ) 
        legend("bottomright",pch=20,col=lcols,legend=levels(torts$Location_ID),cex=0.5)

        layout(matrix(1:4,nrow=2))
        plot( torts$coverage[row(pimat)[ut]], pimat[ut], cex=0.5, pch=20, main='coverage 1' )
        plot( torts$coverage[col(pimat)[ut]], pimat[ut], cex=0.5, pch=20, main='coverage 2' )
        plot( torts$coverage[col(pimat)[ut]]+torts$coverage[row(pimat)[ut]], pimat[ut], cex=0.5, pch=20, main='mean coverage' )
        plot( torts$coverage[col(pimat)[ut]]-torts$coverage[row(pimat)[ut]], pimat[ut], cex=0.5, pch=20, main='diff coverage' )


        ht.lms <- lapply( 1:nrow(pimat), function (kk) {
                lm( pimat[kk,-kk] ~ tort.dists[kk,-kk] )
            } )
        ht.coefs <- sapply( ht.lms, coef )

        layout(1)
        plot( tort.dists[ut], pimat[ut], xlab="geographic distance", ylab="pairwise divergence", pch=20, cex=0.5, col=var_col("Adapter11")[ut] )
        invisible( lapply( ht.lms, function (x) abline(coef(x),col=adjustcolor("black",.1)) ) )

        if (diff(range(diag(pimat)))>0) {
            plot(elev,main="Tortoise heterozygosity")
            lines(county_lines)
            ncols <- 16
            cols <- adjustcolor(diverge_hcl(ncols),.7)
            pifac <- cut( diag(pimat), breaks=ncols )
            points(tort.coords.rasterGCS,pch=21,cex=2,col=grey(.2), bg=cols[pifac])
            # text(tort.coords.rasterGCS,labels=gsub("etort.","",torts$EM_Tort_ID))
            legend("bottomright",pch=20,col=cols,legend=levels(pifac))
        }

        plot(elev,main="Slope in isolation by distance relationship")
        lines(county_lines)
        ncols <- 16
        cols <- adjustcolor(diverge_hcl(ncols),.7)
        slfac <- cut( ht.coefs[,2], breaks=ncols )
        points(tort.coords.rasterGCS,pch=21,cex=2,col=grey(.2), bg=cols[slfac])
        # text(tort.coords.rasterGCS,labels=gsub("etort.","",torts$EM_Tort_ID))
        legend("bottomright",pch=20,col=cols,legend=levels(slfac))


    dev.off()

}
