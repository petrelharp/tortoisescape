library(raster)
library(sp)
library(rgeos)
library(ape)
library(jsonlite)
# library(msarg)
devtools::load_all("~/projects/msarg")
# library(landsim)
devtools::load_all("~/projects/patchy-landscapes/landsim")

# load tortoise habitat
full.habitat <- raster("../visualization/nussear_masked.grd")
habitat <- aggregate( full.habitat, fact=32 )

# population setup
pop <- population( habitat, genotypes='a', accessible=rep(TRUE,length(values(habitat))), 
                  habitable=( !is.na(values(habitat)) & values(habitat)>0 ) )
pop.xy <- xyFromCell(habitat,cell=which(pop$accessible))

# These are the two pairs of twins/repeated samples:
#   etort-143 etort-297
#   etort-156 etort-296
# We will omit these.
rellies <- c("etort-296","etort-297")

# actual sampling locations
sample.coords <- read.csv("../tort_272_info/long-lat.csv",header=TRUE)
sample.coords <- subset(sample.coords, ! (sample.coords$id %in% rellies ) )
sample.points <- spTransform( 
        SpatialPoints( sample.coords[,1:2], 
                      proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") ), 
             CRS(proj4string(pop$habitat)) 
    )
row.names(sample.points) <- sample.ids <- as.character(sample.coords$id)
sample.cells <- cellFromXY( pop$habitat, sample.points )
# the "sample config" to pass to msarg: columns are row, column, layer, #.samples
sample.config <- sort_sample_config( cbind(
                       arrayInd( unique(sample.cells), .dim=as.integer(dim(habitat)[c(2,1,3)]) ),
                       n=table(sample.cells) ) )
# keep track of which sample goes to which simulated sample:
#  sample.points[k] corresponds to the sample.order[k]-th simulated indiv
sample.rowcol <- arrayInd( sample.cells, .dim=as.integer(dim(habitat)[c(2,1,3)]) )
sample.order <- order(sample.rowcol[, 3], sample.rowcol[, 2], sample.rowcol[, 1])
#  and order(sample.order) is the inverse permutation
#  so this is the permutation to re-order the ms output
msord <- order(sample.order)

# observed pairwise divergence
dist.df <- read.csv("../tort_272_info/all_angsd_snps.pwp.csv",header=TRUE,stringsAsFactors=FALSE)
dist.df <- subset( dist.df, ! ( ( etort1 %in% rellies ) | ( etort2 %in% rellies ) ) )
dist.df$etort1 <- factor(dist.df$etort1,levels=sample.ids)
dist.df$etort2 <- factor(dist.df$etort2,levels=sample.ids)
# from Evan 2/10/15
dist.df$years <- dist.df$pi  / 2.064406e-8 / 2
# referenced in www.sciencedirect.com/science/article/pii/S0006320713003443
dist.df$generations <- dist.df$years / 25

# make things necessary to make nice plots of the results
geog.df <- read.csv("../tort_272_info/geog_distance.csv",header=TRUE,stringsAsFactors=FALSE)
geog.df <- subset( geog.df, ! ( ( etort1 %in% rellies ) | ( etort2 %in% rellies ) ) )
geog.df$etort1 <- factor(geog.df$etort1,levels=sample.ids)
geog.df$etort2 <- factor(geog.df$etort2,levels=sample.ids)
dist.df <- merge(dist.df,geog.df)

# use to convert things indexed by location to things indexed by sample
l2s <- rep(1:nrow(sample.config),sample.config[,"n"])
# map tree tips to sample.config
tfn <- tip_order_fn(sample.config)

load("../visualization/north-polygon.RData")  # provides a SpatialPolygon, 'norths'
north.cells <- as.vector(gContains( norths, sample.points, byid=TRUE ))
north.cols <- ifelse( north.cells, "blue", "purple" )
dist.df$col <- ifelse( xor(north.cells[dist.df[,1]],north.cells[dist.df[,2]]), "red",
                      ifelse( north.cells[dist.df[,1]], "blue", "purple" ) )

# of the above, we will need pop, pop.xy and sample.config only.

# helper functions
source("grid-refugia-fns.R")

###
# INITIAL PARAMETERS

init.params <- list(
        pop.density = 0.005/1e4,   # density in indivs/m^2 : one tortoise per 20 ha in perfect habitat; but one tenth this looks better
        sigma = 0.5e3,            # m/gen
        refugia.coords = cbind( x=c(727229,  626639), 
                                y=c(3803747,3957436) ),   # in meters (UTM)
        refugia.radii = 5e4,      # in meters
        refugia.time = 50e4/25,   # time refugia began, in generations
        expansion.time = 20e4/25, # time refugia ended, in generations
        expansion.speed = 0.4,    # in m/gen (not really??)
        expansion.width = 100e3   # in meters, roughly
)

### Other (unchanging) parameters
ntrees <- 100    # Number of trees to simulate.

run.id <- floor(1e6*runif(1))
basedir <- sprintf("run_%06d",run.id)
dir.create( basedir )
cat("Writing out to", basedir, "----\n")

cat( toJSON(init.params, pretty=TRUE), file=file.path(basedir,"params.json") )
write.csv( dist.df, file=file.path(basedir,"dist-df.csv"), row.names=FALSE )

# explore parameter space
# iter.num <- 1
for (iter.num in 1:100) {

    params <- lapply( init.params, function (x) {
                x * 0.5+runif(length(x))
        } )
    params$refugia.coords <- sampleRandom(habitat,size=2,xy=TRUE)[,1:2]

    # (set and) save out the random seed for reproducibility/debugging
    new.seed <- as.integer(runif(1)*2e9)
    set.seed(new.seed)

    outdir <- file.path(basedir, sprintf("iter_%04d",iter.num) )
    dir.create( outdir )
    cat(" ... working on", outdir, "\n")

    cat( toJSON(c(params,list(seed=new.seed)), pretty=TRUE), file=file.path(outdir,"params.json") )

    dem <- do.call( model_setup, c( list(pop=pop), params ) )

    mean.dist <- sim_data( dem, sample.config, outdir, ntrees )
    mean.dist <-  tryCatch(  
                           sim_data( dem, sample.config, outdir, ntrees ), 
                           error=function (e) { 
                               cat("Whoops: skipping this one:\n"); 
                               cat("  ", e$message,"\n"); NULL } )
    if (is.null(mean.dist)) { next }  # skip this one
    rownames(mean.dist) <- colnames(mean.dist) <- sample.ids[sample.order]

    # in the same order as dist.df:
    sim.dist <- mean.dist[msord,msord][ cbind(
                    match(dist.df[,1],sample.ids),
                    match(dist.df[,2],sample.ids) ) ]

    write( sim.dist, file=file.path(outdir,"sim-distances.csv"), ncolumns=1 )

    # mean-square difference:
    model.score <- mean( (dist.df$generations - sim.dist)^2 )

    write( model.score, file=file.path(outdir,"model.score") )

    ### plotting

    # pcs of the covariance matrix (up to scaling; see McVean)
    covmat <- (rowMeans(mean.dist[msord,msord]) + colMeans(mean.dist[msord,msord]) - mean(mean.dist[msord,msord]) - mean.dist[msord,msord])
    pmat <- diag(nrow(covmat)) - 1/nrow(covmat)
    pcs <- eigen( (pmat %*% covmat %*% t(pmat)) )$vectors[,1:2]
    pc.pal <- rainbow(n=32, start=4/6, end=0)
    pc.cols <- apply(pcs, 2, function (x) { pc.pal[cut(x,length(pc.pal))] } )

    # refugia centers
    refugia.centers <- SpatialPoints(coords=params$refugia.coords,
                            proj4string=CRS(proj4string(pop$habitat)))
    refugia <- gBuffer( refugia.centers, width=params$refugia.radii, byid=TRUE )

    # make the plots
    # pdf(file=file.path(outdir,"trees-and-things.pdf"),width=10,height=5,pointsize=10)
    png(file=file.path(outdir, "trees-and-things-%02d.png"), 
        width=10*144, height=5*144, pointsize=10, res=144)
    layout(t(1:2))

    plot( dist.df$distance/1e3, sim.dist, pch=20, cex=0.5, 
         xlab="geog dist (km)", ylab="mean TMRCA (y)",
         col=adjustcolor(dist.df$col,0.5),
         main="simulated vs distance" )
    plot( dist.df$generations, sim.dist, pch=20, cex=0.5, xlim=range(dist.df$generations[dist.df$etort1!=dist.df$etort2]),
         xlab="observed divergence", ylab="mean TMRCA (y)",
         col=adjustcolor(dist.df$col,0.5),
         main="simulated vs observed" )
    dist.lm <- lm( sim.dist ~ dist.df$generations, subset= (dist.df$etort1 != dist.df$etort2) )
    abline(coef(dist.lm))

    # comparisons to individual tortoises
    layout(t(1:3))
    for (tort in paste0("etort-",c(243,35,100,262,218,52,274,90,36))) {
        ut <- with(dist.df,etort1==tort|etort2==tort)
        other <- ifelse( dist.df$etort1[ut]==tort, dist.df$etort2[ut], dist.df$etort1[ut] )
        plot( habitat, main=tort )
        plot( refugia, col=adjustcolor('black',0.25), add=TRUE )
        points( sample.points, pch=ifelse(sample.ids==tort,8,20), cex=2, col=pc.cols[,1] )
        plot( dist.df$distance/1e3, dist.df$generations, pch=20, cex=0.5, ylim=range(dist.df$generations[dist.df$etort1!=dist.df$etort2]),
             xlab="geog dist (km)", ylab="mean TMRCA (y)", main="observed" )
        points( dist.df$distance[ut]/1e3, dist.df$generations[ut], pch=20, cex=2,
             col=pc.cols[other,1] )
        plot( dist.df$distance/1e3, sim.dist, pch=20, cex=0.5,
             xlab="geog dist (km)", ylab="mean TMRCA (y)", main="simulated" )
        points( dist.df$distance[ut]/1e3, sim.dist[ut], pch=20, cex=2,
             col=pc.cols[other,1] )
    }

    layout(t(1:2))
    for (k in 1:ncol(pcs)) {
        plot( habitat, main=sprintf("PC %d",k) )
        plot( refugia, col=adjustcolor('black',0.25), add=TRUE )
        points( sample.points, pch=20, cex=2, col=pc.cols[,k] )
        # if (interactive() && is.null(locator(1))) { break }
    }

    tree.output <- trees_from_ms( file.path(outdir,"msoutput.txt") )
    plot.ntrees <- 4

    for (tree in tree.output[1:plot.ntrees]) {
        plot( habitat, )
        plot( refugia, col=adjustcolor('black',0.25), add=TRUE )
        points( sample.points, pch=20, cex=2, col=pc.cols[,1] )
        # plot_sample_config( dem, sample.config, add=TRUE, xy=pop.xy, col=north.cols )
        pp <- plot.phylo( tree, show.tip.label=FALSE )  # tip.color=pc.cols[,k] )  # tip.color=north.cols, cex=0.2 )
        axisPhylo(1)
        ab <- abline_phylo(v=dem@t[length(dem@t)], lty=2, col='grey')
        points( rep(1.05*ab[1],nrow(pcs)), 1:nrow(pcs), pch=20, col=pc.cols[,1][sample.order][order(tfn(tree))] )
        # if (interactive() && is.null(locator(1))) { break }
    }

    dev.off()
}
