#!/usr/bin/Rscript

require(raster)

sample.loc.obj <- load("geog_coords.RData")
assign("sample.locs",get(sample.loc.obj))
dem <- raster("../visualization/dem_30.gri")

pcs <- read.csv("pcs.csv")

pi <- read.csv("all_angsd_snps.pwp.csv")
# from Evan 1/31/15
pi$pi <- pi$pi * ( 56575857 / 1898838430 ) / 2.01e-9

geog <- read.csv("geog_distance.csv")
geog$distance <- geog$distance / 1000 # now in km

x <- merge( geog, pi, by=c("etort1","etort2") )
x$dpc <- abs(pcs$PC1[match(x$etort1,pcs$etort)] - pcs$PC1[match(x$etort2,pcs$etort)])
x <- subset(x,etort1!=etort2 & (! etort1%in% c("etort-296","etort-297") ) & (! etort2%in% c("etort-296","etort-297") ) )

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

x$group1 <- pop.df$group[match(x$etort1,pop.df$etort)]
x$group2 <- pop.df$group[match(x$etort2,pop.df$etort)]

png(file="ibd-by-pops.png",width=10*288,height=10*288,pointsize=10,res=288)
#pdf(file="ibd-by-pops.pdf",width=10,height=10,pointsize=10)
layout(matrix(1:16,nrow=4))
par(mar=c(3,3,3,1)+.1,mgp=c(2.2,1,0))
plot(dem,legend=FALSE,xaxt='n',yaxt='n')
for (k in seq_along(pops)) { points(pops[[k]],pch=20,col=k,cex=2) }
# legend("topleft",pch=20,cex=2,legend=names(pops),col=seq_along(pops))
# plot( pi ~ distance, data=x, xlab='distance (km)', ylab='divergence', main='all',
#    pch=20, cex=0.25, ylim=range(x$pi), xlim=c(0,max(x$distance)) )
for (j in 1:length(pops)) {
    for (k in j:length(pops)) {
        usethese <- with(x, ( ( group1==names(pops)[j] ) & ( group2==names(pops)[k] ) ) 
            | ( ( group1==names(pops)[k] ) & ( group2==names(pops)[j] ) ) )
        plot( pi ~ distance, data=x, xlab='distance (km)', ylab='divergence', 
            col=adjustcolor("slategrey",.1),
            main=paste(names(pops)[c(j,k)],collapse='-'), 
            pch=20, cex=0.25, ylim=range(x$pi), xlim=c(0,max(x$distance)) )
        points( pi ~ distance, data=subset(x,usethese), xlab='distance (km)', ylab='divergence', 
            pch=20, cex=0.25 )
        points( 3e2+c(-.2e2,0.2e2), c(2.2e6,2.2e6), pch=20, cex=4, col=c(j,k) )
        this.lm <- lm( pi ~ distance, data=subset(x,usethese) )
        abline( coef(this.lm) )
        text( 3e2, 2.3e6, labels=sprintf("y = %2.0f x + %2.0f", coef(this.lm)[2], coef(this.lm)[1]) )
    }
}
dev.off()



png(file="ibd-by-some-pops.png",width=6.5*288,height=2.4*288,pointsize=10,res=288)
layout(t(1:3))
par(mar=c(3,3,3,1)+.1,mgp=c(2.2,1,0))
plot(dem,legend=FALSE,xaxt='n',yaxt='n')
south.pops <- c("west","south","east")
south.col <- "purple"
north.col <- "blue"
ns.cols <- c(north.col,south.col)[1+(names(pops) %in% south.pops)]
for (k in seq_along(pops)) { points(pops[[k]],pch=20,col=ns.cols[k],cex=2) }
both.south <- with(x, ( ( group1%in%south.pops ) & ( group2%in%south.pops ) ) )
both.north <- with(x, ( (! group1%in%south.pops ) & (! group2%in%south.pops ) ) )
# within
    plot( pi ~ distance, data=x, xlab='distance (km)', ylab='mean divergence (years)', 
        col=adjustcolor("slategrey",.1),
        main="within groups",
        pch=20, cex=0.25, ylim=range(x$pi), xlim=c(0,max(x$distance)) )
    points( pi ~ distance, data=subset(x,both.north), xlab='distance (km)', ylab='divergence', pch=20, cex=0.25, col=adjustcolor(north.col,.25) )
    points( pi ~ distance, data=subset(x,both.south), xlab='distance (km)', ylab='divergence', pch=20, cex=0.25, col=adjustcolor(south.col,.25) )
    this.lm <- lm( pi ~ distance, data=subset(x,both.south|both.north) )
    abline( coef(this.lm), col='black' )
    # south.lm <- lm( pi ~ distance, data=subset(x,both.south) )
    # abline( coef(south.lm), col=south.col )
    # north.lm <- lm( pi ~ distance, data=subset(x,both.north) )
    # abline( coef(north.lm), col=north.col )
# between
    plot( pi ~ distance, data=x, xlab='distance (km)', ylab='mean divergence (years)', 
        col=adjustcolor("slategrey",.1),
        main="between groups",
        pch=20, cex=0.25, ylim=range(x$pi), xlim=c(0,max(x$distance)) )
    points( pi ~ distance, data=subset(x,!(both.south|both.north)), xlab='distance (km)', ylab='divergence', pch=20, cex=0.25, col="black" )
    between.lm <- lm( pi ~ distance, data=subset(x,!(both.south|both.north)) )
    abline( coef(between.lm), col='black' )
dev.off()
