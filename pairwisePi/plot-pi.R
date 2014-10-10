#!/usr/bin/Rscript
require(rgdal)
require(maps)
require(maptools)
require(raster)

nind <- 180
read.pimat <- function (filename,tri) {
    pimat.vals <- scan(filename)
    diag <- ( length(pimat.vals)!=choose(nind,2) )
    pimat <- numeric(nind^2); dim(pimat) <- c(nind,nind)
    pimat[tri(pimat,diag=diag)] <- pimat.vals
    pimat[row(pimat)!=col(pimat)] <- (pimat+t(pimat))[row(pimat)!=col(pimat)]
    return( pimat )
}

pimats <- list( 
        robust = read.pimat("alleleCounts_1millionloci.pwp",tri=lower.tri), 
        angsd = read.pimat("350000.pwp",tri=upper.tri), 
        unbiased10 = read.pimat("unbiased/10millionALLsites.pwp",tri=upper.tri), 
        read2 = read.pimat("unbiased/1millionALLsites_exactly2reads.pwp",tri=upper.tri), 
        read3 = read.pimat("unbiased/1millionALLsites_exactly3reads.pwp",tri=upper.tri)
    )

sapply( pimats, mean )

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
dists <- cbind( dists, do.call( cbind, lapply( pimats, function (x) {
                x[ cbind( dists$etort1, dists$etort2 ) ] } ) ) )


# compare the methods
ut <- upper.tri(tort.dists,diag=FALSE) 
pdf(file="pi-methods-comparison.pdf",width=12,height=12,pointsize=10)
lcols <- adjustcolor(rainbow(nlevels(torts$Location_ID)),0.75)
pairs( dists[,-(1:2)], pch=20, cex=0.5, 
        upper.panel=function (x,y,...) points(x,y,col=lcols[torts$Location_ID][row(tort.dists)[ut]],...),
        lower.panel=function (x,y,...) points(x,y,col=lcols[torts$Location_ID][row(tort.dists)[ut]],...)
        )
dev.off()

###
# make plots for each method
for (meth in names(pimats)) {

    pimat <- pimats[[meth]]

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

###
# look at wierdos

lowcoverage <- c("etort-106", "etort-105", "etort-103", "etort-109", "etort-102", "etort-107", "etort-108","etort-104")
highcoverage <- c("etort-61", "etort-127", "etort-151", "etort-141", "etort-31", "etort-16", "etort-51", "etort-110")

pdf(file="comparison-by-coverage.pdf",width=10,height=5,pointsize=10)
layout(t(1:2))

with( dists, plot( unbiased10, read3, pch=20, cex=0.5 ) )
with( subset(dists, ( etort1 %in% lowcoverage ) | ( etort1 %in% lowcoverage ) ), 
    points( unbiased10, read3, col=rainbow(20)[1+ifelse(etort1%in%lowcoverage, match(etort1,lowcoverage),match(etort1,lowcoverage))] ) )
with( subset(dists, ( etort1 %in% highcoverage ) | ( etort1 %in% highcoverage ) ), 
    points( unbiased10, read3, pch=20, col=rainbow(20)[11+ifelse(etort1%in%highcoverage, match(etort1,highcoverage),match(etort1,highcoverage))] ) )
abline(0,1)
legend("bottomright",col=rainbow(20)[c(1:8,11:20)], legend=c(paste("low:",lowcoverage),paste("high:",highcoverage)), pch=c(rep(1,length(lowcoverage)),rep(20,length(highcoverage))) )

with( dists, plot( read2, read3, pch=20, cex=0.5 ) )
with( subset(dists, ( etort1 %in% lowcoverage ) | ( etort1 %in% lowcoverage ) ), 
    points( read2, read3, col=rainbow(20)[1+ifelse(etort1%in%lowcoverage, match(etort1,lowcoverage),match(etort1,lowcoverage))] ) )
with( subset(dists, ( etort1 %in% highcoverage ) | ( etort1 %in% highcoverage ) ), 
    points( read2, read3, pch=20, col=rainbow(20)[11+ifelse(etort1%in%highcoverage, match(etort1,highcoverage),match(etort1,highcoverage))] ) )
abline(0,1)
legend("bottomright",col=rainbow(20)[c(1:8,11:20)], legend=c(paste("low:",lowcoverage),paste("high:",highcoverage)), pch=c(rep(1,length(lowcoverage)),rep(20,length(highcoverage))) )

dev.off()


###
# how is 'unbiased10' without low coverage indivs?

usethese <- with( dists, ! ( ( etort1 %in% lowcoverage ) | ( etort2 %in% lowcoverage ) ) )
tcols <- rainbow(240)
lcols <- adjustcolor(rainbow(nlevels(torts$Location_ID)),0.75)

# these ones still weird:
weirdos <- c("etort-48", "etort-46", "etort-50", "etort-57", "etort-72", "etort-78", "etort-68", "etort-71", "etort-97", "etort-59", "etort-65")

pdf(file="coverage-weirdos.pdf",width=12,height=4)

layout(t(1:2))
plot( torts$coverage, rowMeans(pimats[["unbiased10"]]), type='n', ylab='mean divergence (unbiased10)' )
text( torts$coverage, rowMeans(pimats[['unbiased10']]), labels=torts$EM_Tort_ID, col=1+torts$EM_Tort_ID%in%lowcoverage+2*torts$EM_Tort_ID%in%weirdos  )
plot(elev)
points( tort.coords.rasterGCS, col=1+torts$EM_Tort_ID%in%lowcoverage+2*torts$EM_Tort_ID%in%weirdos, pch=20  )

dev.off()

