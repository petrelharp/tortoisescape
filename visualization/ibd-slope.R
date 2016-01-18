#!/usr/bin/Rscript

require(raster)

sample.loc.obj <- load("geog_coords.RData")
assign("sample.locs",get(sample.loc.obj))
dem <- raster("../geolayers/expanded/16x/dem_30.gri")

pcs <- read.csv("pcs.csv")

pi <- read.csv("all_angsd_snps.pwp.csv",stringsAsFactors=FALSE)
# from Evan 1/31/15
pi$pi <- pi$raw_pi * ( 56575857 / 1898838430 )
pi$generations <- pi$pi  / 2.064406e-8 / 2
# referenced in www.sciencedirect.com/science/article/pii/S0006320713003443
pi$years <- pi$generations * 25

sample.ids <- unique( c(pi$etort1, pi$etort2) )
sample.ids <- sample.ids[ order(as.numeric(gsub(" .*","",gsub("^[^0-9]*","",sample.ids)))) ]
pi$etort1 <- factor(pi$etort1,levels=sample.ids)
pi$etort2 <- factor(pi$etort2,levels=sample.ids)

geog <- read.csv("geog_distance.csv")
geog$distance <- geog$distance / 1000 # now in km

torts <- read.csv("sample_metadata.csv",stringsAsFactors=FALSE)
torts$etort <- factor( torts$EM_Tort_ID , levels=sample.ids )

dists <- merge( geog, pi, by=c("etort1","etort2") )
dists$dpc <- abs(pcs$PC1[match(dists$etort1,pcs$etort)] - pcs$PC1[match(dists$etort2,pcs$etort)])
dists <- subset(dists,etort1!=etort2 & (! etort1%in% c("etort-296","etort-297") ) & (! etort2%in% c("etort-296","etort-297") ) )

coords <- coordinates(sample.locs)
pops <- list(
        west = sample.locs[ 
            ( coords[,1] < -18e5 ) & 
            ( coords[,2] > -3e5 ) & 
            ( coords[,2] < -1.5e5 ) 
        ],
        south = sample.locs[
            ( coords[,2] < -3e5 ) &
            ( coords[,1]+coords[,2] < -21e5 ) 
        ],
        east = sample.locs[
            ( coords[,1]+coords[,2] > -21e5 ) &
            ( coords[,1]-coords[,2] > -15.1e5 ) &
            ( pcs$PC1 < 0 )  # ( coords[,2] < -2e5 )
        ],
        north_east = sample.locs[
            ( coords[,2]+coords[,1]/5 > -4.7e5 )
        ]
    )
pops$north <- sample.locs[ setdiff( row.names(sample.locs), unlist(lapply(pops,row.names)) ) ]

pop.df <- data.frame( 
        etort=unlist(lapply(pops,row.names)),
        group=names(pops)[rep(seq_along(pops),unlist(lapply(pops,length)))],
    stringsAsFactors=FALSE )

dists$group1 <- pop.df$group[match(dists$etort1,pop.df$etort)]
dists$group2 <- pop.df$group[match(dists$etort2,pop.df$etort)]

png(file="ibd-by-pops.png",width=10*288,height=10*288,pointsize=10,res=288)
#pdf(file="ibd-by-pops.pdf",width=10,height=10,pointsize=10)
pifac <- 1000  # in kb
layout(matrix(1:16,nrow=4))
par(mar=c(3,3,3,1)+.1,mgp=c(2.2,1,0))
plot(dem,legend=FALSE,xaxt='n',yaxt='n')
for (k in seq_along(pops)) { points(pops[[k]],pch=20,col=k,cex=2) }
posfn <- function (eps) { pifac*((1-eps)*min(dists$pi) + eps*max(dists$pi)) }
for (j in 1:length(pops)) {
    for (k in j:length(pops)) {
        usethese <- with(dists, ( ( group1==names(pops)[j] ) & ( group2==names(pops)[k] ) ) 
            | ( ( group1==names(pops)[k] ) & ( group2==names(pops)[j] ) ) )
        plot( pifac*pi ~ distance, data=dists, xlab='distance (km)', ylab='divergence', 
            col=adjustcolor("slategrey",.1),
            main=paste(names(pops)[c(j,k)],collapse='-'), 
            pch=20, cex=0.25, ylim=pifac*range(dists$pi), xlim=c(0,max(dists$distance)) )
        points( pifac*pi ~ distance, data=subset(dists,usethese), xlab='distance (km)', ylab='divergence', 
            pch=20, cex=0.25 )
        points( 3e2+c(-.2e2,0.2e2), posfn(c(.05,.05)), pch=20, cex=4, col=c(j,k) )
        this.lm <- lm( pi ~ distance, data=subset(dists,usethese) )
        abline( pifac*coef(this.lm) )
        text( 3e2, posfn(.15), labels=sprintf("y = %0.4f x + %0.4f", pifac*coef(this.lm)[2], pifac*coef(this.lm)[1]) )
    }
}
dev.off()

require(rgdal)
county.lines.obj <- load("../visualization/county_lines.Robj") 
assign("county.lines", spTransform(get(county.lines.obj),  CRSobj=CRS(proj4string(dem)) ) )
roads <- spTransform( readOGR("../visualization","moj_maj_road"), CRSobj=CRS(proj4string(dem)) )

png(file="ibd-by-some-pops.png",width=6.5*288,height=2.4*288,pointsize=10,res=288)
layout(t(1:3))
par(mar=c(3,3,3,1)+.1,mgp=c(2.2,1,0))
plot(dem,legend=FALSE,xaxt='n',yaxt='n',main='sample locations')
lines(county.lines,col=adjustcolor("slategrey",0.5))
lines(roads, col=adjustcolor("darkslategrey",0.5),lwd=0.5)
south.pops <- c("west","south","east")
south.col <- "purple"
north.col <- "blue"
ns.cols <- c(north.col,south.col)[1+(names(pops) %in% south.pops)]
for (k in seq_along(pops)) { points(pops[[k]],pch=20,col=ns.cols[k],cex=1) }
both.south <- with(dists, ( ( group1%in%south.pops ) & ( group2%in%south.pops ) ) )
both.north <- with(dists, ( (! group1%in%south.pops ) & (! group2%in%south.pops ) ) )
pifac <- 1000  # in kb
# within
posfn <- function (eps) { pifac*( (1-eps)*min(dists$pi) + eps*max(dists$pi) ) }
    plot( pifac*pi ~ distance, data=dists, xlab='distance (km)', ylab='mean divergence (per Kb)', 
        col=adjustcolor("slategrey",.1),
        main="within groups",
        pch=20, cex=0.25, ylim=pifac*range(dists$pi), xlim=c(0,max(dists$distance)) )
    points( pifac*pi ~ distance, data=subset(dists,both.north), pch=20, cex=0.25, col=adjustcolor(north.col,.25) )
    points( pifac*pi ~ distance, data=subset(dists,both.south), pch=20, cex=0.25, col=adjustcolor(south.col,.25) )
    this.lm <- lm( pi ~ distance, data=subset(dists,both.south|both.north) )
    abline( pifac*coef(this.lm), col='black' )
    text( 3e2, posfn(.1), labels=sprintf("y = %0.4f x + %0.4f", pifac*coef(this.lm)[2], pifac*coef(this.lm)[1]) )
# between
    plot( pifac*pi ~ distance, data=dists, xlab='distance (km)', ylab='mean divergence (years)', 
        col=adjustcolor("slategrey",.1),
        main="between groups",
        pch=20, cex=0.25, ylim=pifac*range(dists$pi), xlim=c(0,max(dists$distance)) )
    points( pifac*pi ~ distance, data=subset(dists,!(both.south|both.north)), xlab='distance (km)', ylab='divergence', pch=20, cex=0.25, col="black" )
    this.lm <- lm( pi ~ distance, data=subset(dists,!(both.south|both.north)) )
    abline( pifac*coef(this.lm), col='black' )
    text( 3e2, posfn(.1), labels=sprintf("y = %0.4f x + %0.4f", pifac*coef(this.lm)[2], pifac*coef(this.lm)[1]) )
dev.off()

png(file="ibd-by-some-pops-years.png",width=6.5*288,height=2.4*288,pointsize=10,res=288)
layout(t(1:3))
par(mar=c(3,3,3,1)+.1,mgp=c(2.2,1,0))
plot(dem,legend=FALSE,xaxt='n',yaxt='n')
lines(county.lines,col=adjustcolor("slategrey",0.5))
lines(roads, col=adjustcolor("darkslategrey",0.5),lwd=0.5)
south.pops <- c("west","south","east")
south.col <- "purple"
north.col <- "blue"
ns.cols <- c(north.col,south.col)[1+(names(pops) %in% south.pops)]
for (k in seq_along(pops)) { points(pops[[k]],pch=20,col=ns.cols[k],cex=1) }
both.south <- with(dists, ( ( group1%in%south.pops ) & ( group2%in%south.pops ) ) )
both.north <- with(dists, ( (! group1%in%south.pops ) & (! group2%in%south.pops ) ) )
# within
posfn <- function (eps) { (1-eps)*min(dists$years) + eps*max(dists$years) }
    plot( years ~ distance, data=dists, xlab='distance (km)', ylab='mean divergence (years)', 
        col=adjustcolor("slategrey",.1),
        main="within groups",
        pch=20, cex=0.25, ylim=range(dists$years), xlim=c(0,max(dists$distance)) )
    points( years ~ distance, data=subset(dists,both.north), xlab='distance (km)', ylab='divergence', pch=20, cex=0.25, col=adjustcolor(north.col,.25) )
    points( years ~ distance, data=subset(dists,both.south), xlab='distance (km)', ylab='divergence', pch=20, cex=0.25, col=adjustcolor(south.col,.25) )
    this.lm <- lm( years ~ distance, data=subset(dists,both.south|both.north) )
    abline( coef(this.lm), col='black' )
    text( 3e2, posfn(.1), labels=sprintf("y = %2.0f x + %2.0f", coef(this.lm)[2], coef(this.lm)[1]) )
    # south.lm <- lm( years ~ distance, data=subset(dists,both.south) )
    # abline( coef(south.lm), col=south.col )
    # north.lm <- lm( years ~ distance, data=subset(dists,both.north) )
    # abline( coef(north.lm), col=north.col )
# between
    plot( years ~ distance, data=dists, xlab='distance (km)', ylab='mean divergence (years)', 
        col=adjustcolor("slategrey",.1),
        main="between groups",
        pch=20, cex=0.25, ylim=range(dists$years), xlim=c(0,max(dists$distance)) )
    points( years ~ distance, data=subset(dists,!(both.south|both.north)), xlab='distance (km)', ylab='divergence', pch=20, cex=0.25, col="black" )
    this.lm <- lm( years ~ distance, data=subset(dists,!(both.south|both.north)) )
    abline( coef(this.lm), col='black' )
    text( 3e2, posfn(.1), labels=sprintf("y = %2.0f x + %2.0f", coef(this.lm)[2], coef(this.lm)[1]) )
dev.off()



dups <- list(
        "etort-156"=c("etort-156","etort-296"),
        "etort-143"=c("etort-143","etort-297")
    )

resids <- lapply( dups, function (xx) {
        d1 <- with( subset(pi,etort1==xx[1] | etort2==xx[1]), data.frame(
                etort=factor(sample.ids[ifelse(etort1==xx[1],etort2,etort1)],levels=sample.ids),
                pi1=pi ) )
        d2 <- with( subset(pi,etort1==xx[2] | etort2==xx[2]), data.frame(
                etort=factor(sample.ids[ifelse(etort1==xx[2],etort2,etort1)],levels=sample.ids),
                pi2=pi ) )
        d12 <- merge(d1,d2)
        return( d12 )
    } )

lapply( resids, function (z) { summary( z$pi2-z$pi1 ) } )
lapply( resids, function (z) { sd( z$pi2-z$pi1 ) } )

pdf(file="duplicates.pdf",width=6,height=6,pointsize=10)
layout(matrix(1:4,nrow=2))
lapply( seq_along(resids), function (k) { 
        plot( resids[[k]]$pi1,resids[[k]]$pi2, xlab=dups[[k]][1], ylab=dups[[k]][2] )
        abline(0,1,col='red')
        hist( resids[[k]]$pi2-resids[[k]]$pi1, breaks=40, xlab=paste("residuals,",names(resids)[k]), main="" ) 
        abline(v=0,col='red',lwd=2)
    } )
dev.off()


##
# only within Ivanpah

torts$Location_ID <- factor(torts$Location_ID)
ivanpah.loc.names <- c("Silver State", "ISEGS", "Ivanpah")
ivanpah.torts <- torts$EM_Tort_ID[ (torts$Location_ID %in% ivanpah.loc.names) & ! (torts$EM_Tort_ID %in% sapply(dups,"[",1)) ]
ivanpah.locs <- sample.locs[ match(ivanpah.torts,row.names(sample.locs)) ]

ivanpah.loc.factor <- levels(torts$Location_ID)[ torts$Location_ID[ torts$Location_ID %in% ivanpah.loc.names ] ]
ivanpah.loc.factor <- factor( paste( ivanpah.loc.factor, ifelse( ivanpah.loc.factor=="ISEGS", 
                        ifelse( torts$Easting[torts$Location_ID %in% ivanpah.loc.names] < 645500, "_east", "_west" ),
                    "" ), sep='' ) )

ivanpah.dists <- subset( dists, ( etort1 %in% ivanpah.torts ) & ( etort2 %in% ivanpah.torts ) & ( etort1 != etort2 ) )
ivanpah.dists$loc1 <- ivanpah.loc.factor[ match(ivanpah.dists$etort1,ivanpah.torts) ]
ivanpah.dists$loc2 <- ivanpah.loc.factor[ match(ivanpah.dists$etort2,ivanpah.torts) ]
ivanpah.dists$col <- with(ivanpah.dists, ifelse( loc1==loc2, "green", "black") )

ivanpah.lm <- lm( pi ~ distance, data=ivanpah.dists )
ivanpah.within.lm <- lm( pi ~ distance, data=subset(ivanpah.dists,loc1==loc2) )

infl <- function (x,fac=.5) { as.vector(mean(x) + (1+fac)*(x-mean(x))) }

require(plotrix)

png(file="within-ivanpah-ibd.png",res=200,width=6*200,height=4*200)
layout(t(1:2))
yfac <- 1000
    posfn <- function (eps) { yfac * ( (1-eps)*min(ivanpah.dists$pi) + eps*max(ivanpah.dists$pi) ) }
    opar <- par(mar=c(5,4,3,1)+.1)
    with(ivanpah.dists, plot( distance, yfac*pi, 
			pch=20, cex=0.3, 
            xlab="geographic distance (km)", ylab="genetic distance (per Kb)") )
    title( main="Isolation by Distance", line=1 )
    abline( yfac*coef( ivanpah.lm ), col='red' )
    text( 25, posfn(.1), labels=sprintf("y = %0.4f x + %0.4f", yfac*coef(ivanpah.lm)[2], yfac*coef(ivanpah.lm)[1]) )
    par(mar=c(5,0,3,1)+.1)
    plot(dem, xlim=infl(ivanpah.locs@bbox["coords.x1",],1.2), ylim=infl(ivanpah.locs@bbox["coords.x2",]), xaxt='n', yaxt='n', 
        smallplot=c(0.68,0.71,0.2845,0.8155), 
        bigplot=c(0,0.64,0.255,0.845), 
        legend.lab="elevation (m)" )
    title( main="in the Ivanpah valley", line=1 )
    points( ivanpah.locs, pch=20 ) #, col=torts$Location_ID[torts$Location_ID %in% ivanpah.loc.names] )
    scalebar(1e4,label=c("","10km",""))
dev.off()
