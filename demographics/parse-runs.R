# Get the results from all the iterations
library(jsonlite)

# function to examine output of a single run
get_results <- function (basedir=".") {
    all.iter.dirs <- grep( "run_.*iter_", list.dirs( basedir ), value=TRUE )
    have.score <- file.exists( file.path(all.iter.dirs,"model.score") )
    iter.dirs <- all.iter.dirs[have.score]
    base.dirs <- unique( dirname(all.iter.dirs) )
    base.params.list <- lapply( base.dirs, function (itd) { fromJSON(file.path(itd,"params.json")) } )
    base.params <- do.call( rbind, lapply( base.params.list, function (x) {
                        dim(x$refugia.coords) <- c(1,4)
                        if (is.null(x$hab.fact)) { x$hab.fact <- 32 }
                        as.data.frame(x)
        } ) )
    base.params$dir <- base.dirs
    all.params.list <- lapply( iter.dirs, function (itd) { fromJSON(file.path(itd,"params.json")) } )
    all.params <- do.call( rbind, lapply( all.params.list, function (x) {
                        dim(x$refugia.coords) <- c(1,4)
                        as.data.frame(x)
            } ) )
    all.params$score <- sapply( file.path(iter.dirs,"model.score"), scan, quiet=TRUE )
    all.params$dir <- iter.dirs
    base.match <- match( dirname(all.params$dir), base.params$dir  )
    all.params$ntrees <- base.params$ntrees[base.match]
    all.params$hab.fact <- base.params$hab.fact[base.match]
    all.params <- all.params[ order(all.params$score,decreasing=FALSE), ]
    return(all.params)
}

new.all.params <- get_results()
last.all.params <- if (file.exists("all-results.csv")) { read.csv("all-results.csv", header=TRUE) } else { NULL }
all.params <- merge( new.all.params, last.all.params, all.x=TRUE, all.y=TRUE )
all.params <- all.params[ order(all.params$score,decreasing=FALSE), ]

write.csv(all.params, file="all-results.csv", row.names=FALSE)

# pairs(all.params)

# look at variation between runs with the same parameter

# ntrees=500 for these -- looks good.
same <- get_results("run_all_same")
sd(same$score) #= 163433
sd(same$score)/mean(same$score)  # = .037


# save out the good params to start new things from
goodones <- ( all.params$score < 25e4 )
write.csv( subset(all.params, goodones), file="good-results-lt-25e4.csv", row.names=FALSE )

# plot everything so far

# refugia locations
full.habitat <- raster("../visualization/nussear_masked.grd")
r1 <- SpatialPoints( cbind( jitter(all.params$refugia.coords.1), jitter(all.params$refugia.coords.3) ), proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") ) 
r2 <- SpatialPoints( cbind( jitter(all.params$refugia.coords.2), jitter(all.params$refugia.coords.4) ), proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") ) 

# make plots
## holy crap!  there's definately a good region!
param.names <- c("pop.density", "sigma", "refugia.coords.1", "refugia.coords.2", "refugia.coords.3", "refugia.coords.4", "refugia.radii", "refugia.time", "expansion.time", "expansion.speed", "expansion.width", "score")

png( file="all-parameter-pairs.png", width=12*144, height=12*144, res=144, pointsize=10 )
pairs( all.params[c(which(!goodones),which(goodones)),param.names], 
      col=ifelse( goodones, 'red', 'black' )[c(which(!goodones),which(goodones))] )
dev.off()

png( file="all-refugia-pairs.png", width=8*144, height=4*144, res=144, pointsize=10 )
layout( t(1:2) )
plot(full.habitat)
points(r1, col='black', pch=20)
points(r1[goodones], col='red', pch=20)
plot(full.habitat)
points(r2, col='black', pch=20)
points(r2[goodones], col='red', pch=20)

plot( all.params[,c("refugia.coords.1","refugia.coords.3")], 
     col=ifelse( goodones, 'red', 'black' ) )
points( all.params[goodones,c("refugia.coords.1","refugia.coords.3")], 
     col="red", pch=20 )
plot( all.params[,c("refugia.coords.2","refugia.coords.4")], 
     col=ifelse( goodones, 'red', 'black' ) )
points( all.params[goodones,c("refugia.coords.2","refugia.coords.4")], 
     col="red", pch=20 )
dev.off()

png( file="all-parameter-hists.png", width=4*144, height=16*144, res=144, pointsize=10 )
layout(seq_along(param.names))
for (x in param.names) {
    hist( all.params[[x]][!goodones], breaks=50, xlim=range(all.params[[x]]), freq=FALSE, main=x )
    hist( all.params[[x]][goodones], breaks=50, xlim=range(all.params[[x]]), freq=FALSE, 
         add=TRUE, col=adjustcolor("red",0.5) )
}
dev.off()

