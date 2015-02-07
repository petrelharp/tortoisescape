#!/usr/bin/Rscript

require(raster)

sample.loc.obj <- load("geog_coords.RData")
assign("sample.locs",get(sample.loc.obj))
dem <- raster("../visualization/dem_30.gri")

pcs <- read.csv("pcs.csv")
geog <- read.csv("geog_distance.csv")
pi <- read.csv("all_angsd_snps.pwp.csv")

x <- merge( geog, pi, by=c("etort1","etort2") )
x$dpc <- abs(pcs$PC1[match(x$etort1,pcs$etort)] - pcs$PC1[match(x$etort2,pcs$etort)])
x <- subset(x,etort1!=etort2 & (! etort1%in% c("etort-296","etort-297") ) & (! etort2%in% c("etort-296","etort-297") ) )

coords <- coordinates(sample.locs)
pops <- list(
        west_mojave = sample.locs[ 
            ( coords[,1] < -18e5 ) & 
            ( coords[,2] > -3e5 ) & 
            ( coords[,2] < -1.5e5 ) 
        ],
        south_mojave = sample.locs[
            ( coords[,2] < -3e5 ) &
            ( coords[,1]+coords[,2] < -21e5 ) 
        ],
        east_mojave = sample.locs[
            ( coords[,1]+coords[,2] > -21e5 ) &
            ( coords[,1]-coords[,2] > -15.1e5 ) &
            ( pcs$PC1 < 0 )  # ( coords[,2] < -2e5 )
        ]
    )
pops$north_mojave <- sample.locs[ setdiff( row.names(sample.locs), unlist(lapply(pops,row.names)) ) ]

pop.df <- data.frame( 
        etort=unlist(lapply(pops,row.names)),
        group=names(pops)[rep(seq_along(pops),unlist(lapply(pops,length)))],
    stringsAsFactors=FALSE )

x$group1 <- pop.df$group[match(pop.df$etort,x$etort1)]
x$group2 <- pop.df$group[match(pop.df$etort,x$etort2)]

plot(dem)
for (k in seq_along(pops)) {
    points(pops[[k]],pch=20,col=k,cex=2)
}
legend("topleft",pch=20,cex=2,legend=names(pops),col=seq_along(pops))
