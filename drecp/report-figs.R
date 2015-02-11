#!/usr/bin/Rscript

require(raster)
require(rgdal)

dem <- raster("../geolayers/nussear/habitat-model/dem_30.grd")
nus <- raster("../geolayers/nussear/habitat-model/nussear.grd")
crew <- raster("../geolayers/nussear/habitat-model/crew_NA_1km.gri")
bb <- extent( 3.4e5, 8.7e5, 3.62e6, 4.17e6 )

tort.loc.obj <- load("../tort_272_info/geog_coords.RData")
assign( "sample.locs", spTransform( get(tort.loc.obj), CRSobj=CRS(proj4string(nus)) ) )
ref.pt.obj <- load("../geolayers/reference_points/all_ref_points.RData")
assign( "ref.points", spTransform( get(ref.pt.obj), CRSobj=CRS(proj4string(nus)) ) )
isolated.refs <- c( 2, 16, 44, 50, 136, 142, 153, 174 )

png(width=3*288, height=3*288, res=288, file="reference-points.png", pointsize=10)
#pdf(width=2.5, height=2.5, file="reference-points.pdf", pointsize=10)
par(mar=c(1,1,3,1)+.1)
plot( crop(mask(nus,crew),bb), legend=FALSE,xaxt='n',yaxt='n',main='reference points')
points(ref.points[-isolated.refs], pch=20, cex=0.25)
dev.off()

require(knitr)
chunks <- read.csv("chunk-summaries.csv",header=TRUE)
x <- chunks[order(chunks$mean,decreasing=TRUE),c("chunk","area","mean","max.incr")]
x$mean <- round(x$mean)
x$max.incr <- round(x$max.incr)
colnames(x) <- c("chunk","isolated area","mean isolation","maximum isolation")
x <- x[c(match("all",x[,1]),setdiff(1:nrow(x),match("all",x[,1]))),]

kable( x, row.names=FALSE )


####
# parse different alternatives

alt.names <- c( alt_pref="alt_pref_pda", alt_1="alt_1_pda", alt_2="alt_2_pda", alt_3="alt_3_pda", alt_4="alt_4_pda" )
alt.full.names <- c( alt_pref="preferred", alt_1="alt #1", alt_2="alt #2", alt_3="alt #3", alt_4="alt #4" )
alt.envs <- lapply( alt.names, function (aname) {
        aenv <- new.env()
        load(paste("../inference/habitat-only/nus_gt_three/",aname,"_habitat-only_all_chunks.RData",sep=''),envir=aenv)
        assign( "alt", crop( extend( raster(file.path("../geolayers/alternatives",paste(aname,".tif",sep=''))), nalayer ), nalayer ), aenv )
        return(aenv)
    } )

sym <- function (x) { (x+t(x))/2 }
get.nearby <- function (sub.hts,thresh=.4) {
    nearby <- ( sym(sub.hts) < apply(sym(sub.hts),1,quantile,thresh) )  # average over nearest (thresh) locations
    return(nearby)
}

nalayer <- raster("../geolayers/nussear/habitat-model/mask_crew_dem_2K_sea_habitat.grd")
if (!file.exists("ref_point_area.csv")) {
    dpts <- apply( sapply( 1:length(ref.points), function (k) { values( distanceFromPoints(nalayer,ref.points[k]) ) } ), 1, which.min )
    point.area <- tapply( values(nalayer), dpts, sum, na.rm=TRUE )
    write.csv( cbind( ref=seq_along(ref.points), area=point.area ), file="ref_point_area.csv", row.names=FALSE )
} else {
    point.area.df <- read.csv("ref_point_area.csv")
    point.area <- point.area.df$area
}

### table

isolation.thresh <- 15000
nearby.thresh <- 0.4
total.area <- sum(!is.na(values(nalayer)))
alt.stats <- sapply(alt.envs, function (x) { with(x, {
            diffmat <- sym(sub.alt.hts-sub.hts)
            nearby <- get.nearby( sub.hts, nearby.thresh )
            blocked <- attr(sub.alt.hts,"blocked")
            nearby.means <- rowSums( diffmat*nearby ) / rowSums(nearby)
            data.frame(
                "habitat area removed (km2)" = sum( !is.na( values(nalayer)[ unique(c(blocked,which(values(alt)%in%this.chunk))) ] ) ),
                "percent habitat removed" = sum( !is.na( values(nalayer)[ unique(c(blocked,which(values(alt)%in%this.chunk))) ] ) )/total.area,
                "significantly isolated area (km2)" = sum( point.area[nearby.means>isolation.thresh] ),
                "percent significantly isolated area" = sum( point.area[nearby.means>isolation.thresh] )/total.area,
                "mean increase on significantly isolated area" = mean(nearby.means[nearby.means>isolation.thresh]),
                "isolated area (km2)" = sum( point.area[nearby.means>0] ),
                "percent isolated area" = sum( point.area[nearby.means>0] )/total.area,
                "mean increase on isolated area" = mean(nearby.means[nearby.means>0]),
                "max increase" = max( diffmat[nearby], na.rm=TRUE ),
                "max decrease" = min( diffmat[nearby], na.rm=TRUE ),
                "mean relative increase" = mean( ( rowSums( diffmat*nearby/sym(1+sub.hts) )/rowSums(nearby) )[nearby.means>0]),
                "mean relative increase on isolated area" = mean( ( rowSums( diffmat*nearby/sym(1+sub.hts) )/rowSums(nearby) )[nearby.means>isolation.thresh]),
                check.names=FALSE)
        } ) 
    } )
write.csv(alt.stats,file="alternative-summaries.csv")

alt.stats <- read.csv("alternative-summaries.csv",header=TRUE,row.names=1)

alts <- data.frame(
        "habitat removed (km2)"= floor(unlist(alt.stats["habitat area removed (km2)",])),
        "(%)"=sprintf("(%0.1f%%)",
                (100*unlist(alt.stats["percent habitat removed",]))
            ),
        "isolated (km2)"= floor(unlist(alt.stats["isolated area (km2)",])),
        "(%)"=sprintf("(%0.1f%%)",
                (100*unlist(alt.stats["percent isolated area",]))
            ),
        "isolation (years)"= floor(unlist(alt.stats["mean increase on isolated area",])),
        "strongly isolated (km2)"= floor(unlist(alt.stats["significantly isolated area (km2)",])),
        "(%)"=sprintf("(%0.1f%%)",
                (100*unlist(alt.stats["percent significantly isolated area",]))
            ),
        "strong isolation (years)"= floor(unlist(alt.stats["mean increase on isolated area",])),
        "relative isolation "=sprintf("%i%%",
                floor(100*unlist(alt.stats["mean relative increase on isolated area",]))
            ),
        check.names=FALSE
    )
rownames(alts) <- c("preferred",paste("alternative",1:4))

kable(alts)


### figures
require(colorspace)
require(fields)
diff.cols <- diverge_hcl(128, h=c(225,0), c=100, l=c(60,95), power=0.65)
interp.values <- function (x,refs) {
    # interpolate values seen at ref.points[refs] to everywhere else
    stopifnot( ( is.integer(refs) && length(x)==length(refs) ) || ( is.logical(refs) && length(x)==sum(refs) ) )
    tps <- fastTps( coordinates(ref.points)[refs,], x, lambda=1e-8, theta=2e5 )
    return( mask(interpolate(nalayer,tps),nalayer) )
}

tdiffs <- sapply( alt.envs, function (aenv) {
        with(aenv, {
                nearby <- get.nearby(sub.hts)
                x <- rep(NA,length(alt.good.refs))
                x[alt.good.refs] <- ( rowSums( (sym(sub.alt.hts-sub.hts) * nearby), na.rm=TRUE ) / rowSums(nearby) )
                return(x)
            } )
    } )
zlims <- c(-1,1)*max(abs(tdiffs),na.rm=TRUE)

reldiffs <- sapply( alt.envs, function (aenv) {
        with(aenv, {
                nearby <- get.nearby(sub.hts)
                x <- rep(NA,length(alt.good.refs))
                x[alt.good.refs] <- ( rowSums( ( (sym(sub.alt.hts-sub.hts) * nearby) / (sym(sub.hts)) ), na.rm=TRUE ) / rowSums(nearby) )
                return(x)
            } )
    } )
ratio.zlims <- c(-1,1)*max(abs(reldiffs),na.rm=TRUE)

# figures for each alternative
for (k in seq_along(alt.names)) {
    png(file=paste(alt.names[k],"-diffs.png",sep=''),width=6*288,height=3*288,pointsize=10,res=288)
    layout(t(1:2))
    par(mar=c(1,1,2,1)+.1)
    with( alt.envs[[k]], {
            layout(t(1:2))
            x <- tdiffs[,k]
            plot( interp.values(x[alt.good.refs],alt.good.refs), main="mean difference nearby", col=diff.cols, zlim=zlims, 
                legend.width=2, xaxt='n', yaxt='n', legend.mar=8.1 )
            plot( alt, col="slategrey", add=TRUE, legend=FALSE )
            mtext( alt.full.names[k], side=3, line=-1 )
            y <- reldiffs[,k]
            plot( interp.values(y[alt.good.refs],alt.good.refs), main="relative difference nearby", col=diff.cols, zlim=ratio.zlims, 
                legend.width=2, xaxt='n', yaxt='n', legend.mar=5.1 )
            plot( alt, col="lightslategrey", add=TRUE, legend=FALSE )
        } )
    dev.off()
}
