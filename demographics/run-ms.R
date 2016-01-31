source("msarg/msarg.R",chdir=TRUE)
library(raster)
library(rgeos)
library(sp)
library(ape)

# plotting stuff
habitat <- raster("../visualization/nussear_masked.grd")
load("../visualization/north-polygon.RData")  # provides a SpatialPolygon, 'norths'

demogs <- list( 
               refugia.ms="refugia-demography-msarg.RData",
               barrier.ms="barrier-demography-msarg.RData"
              )

ntrees <- 1000

# # first the barrier
# load("barrier-demography-msarg.RData")  # provides barrier.demog and sample.config
# outdir <- "barrier.ms"

for (kdem in seq_along(demogs)) {
    outdir <- names(demogs)[kdem]
    demog.obj <- setdiff( load(demogs[[kdem]]), c("sample.config", "pop.xy") )
    dem <- get(demog.obj)

    # sample locations
    sample.xy <- pop.xy[sample_config_index( sample.config, dem ),]
    sample.points <- SpatialPoints(sample.xy,proj4string=CRS(proj4string(habitat)))
    north.cells <- as.vector(gContains( norths, sample.points, byid=TRUE ))
    geog.dist <- sqrt( outer(sample.xy[,1],sample.xy[,1],"-")^2 + outer(sample.xy[,2],sample.xy[,2],"-")^2 )
    # use to convert things indexed by location (as those above) to things indexed by sample
    l2s <- rep(1:nrow(sample.config),sample.config[,"n"])

    sample.cols <- ifelse( north.cells, "blue", "purple" )
    comparison.cols <- outer( north.cells, north.cells, function (x,y) { ifelse( xor(x,y), "red", ifelse(x, "blue", "purple") ) } )

    # actually run ms
    run_ms( dem, nsamp=sample.config, trees=TRUE, outdir=outdir, nreps=ntrees )

    tree.output <- trees_from_ms( file.path(outdir,"msoutput.txt") )
    tree.dists <- tree_dists( tree.output, sample.config=sample.config )
    mean.dist <- tree.dists[[1]]
    for (k in seq_along(tree.dists)[-1]) { mean.dist <- mean.dist+tree.dists[[k]] }
    mean.dist <- mean.dist/length(tree.dists)

    # write out things necessary to make nice plots of the results
    sample.ids <- paste0("sample-",1:sum(sample.config[,"n"]))
    ut <- upper.tri(mean.dist,diag=FALSE)
    dist.table <- data.frame( id1=sample.ids[row(mean.dist)[ut]], id2=sample.ids[col(mean.dist)[ut]], pi=mean.dist[ut] )
    write.csv( dist.table, file=file.path(outdir,"mean-divergence.csv"), row.names=FALSE )
    geog.table <- data.frame( id1=sample.ids[row(mean.dist)[ut]], id2=sample.ids[col(mean.dist)[ut]], distance=geog.dist[l2s,l2s][ut] )
    write.csv( geog.table, file=file.path(outdir,"geog-distance.csv"), row.names=FALSE )
    sample.spatialpoints <- sample.points[l2s]
    row.names(sample.spatialpoints) <- sample.ids
    save( sample.spatialpoints, file=file.path(outdir,"geog_coords.RData"))
    fake.pcs <- data.frame( id=sample.ids, pc1=ifelse(north.cells[l2s],+1.0,-1.0) )
    write.csv( fake.pcs, file=file.path(outdir,"pcs.csv"), row.names=FALSE )

    ## look at some plots

    pdf(file=file.path(outdir,"trees-and-things.pdf"),width=10,height=5,pointsize=10)

    plot( as.vector(geog.dist[l2s,l2s])/1e3, as.vector(mean.dist), pch=20, cex=0.5, 
         xlab="geog dist (km)", ylab="mean TMRCA (y)",
         col=adjustcolor(comparison.cols[l2s,l2s],0.5) )

    layout(t(1:2))
    for (tree in tree.output) {
        plot(habitat)
        plot_sample_config( dem, sample.config, add=TRUE, xy=pop.xy, cex=0.1, sample.cols=sample.cols )
        plot.phylo( tree, tip.color=sample.cols[l2s], cex=0.2 )
        ape::axisPhylo(1)
        abline_phylo(v=dem@t, lty=2, lwd=2, col='grey')
        if (interactive() && is.null(locator(1))) { break }
    }

    dev.off()
}
