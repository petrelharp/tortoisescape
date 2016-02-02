library(raster)
library(rgeos)
library(sp)
library(ape)
library(msarg)

# plotting stuff
habitat <- raster("../visualization/nussear_masked.grd")
load("../visualization/north-polygon.RData")  # provides a SpatialPolygon, 'norths'

demogs <- list( 
               refugia.ms="refugia-demography-msarg.RData",
               barrier.ms="barrier-demography-msarg.RData"
              )

ntrees <- 500
plot.ntrees <- min(20,ntrees)

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
    # map tree tips to sample.config
    tfn <- tip_order_fn(sample.config)

    sample.cols <- ifelse( north.cells, "blue", "purple" )
    comparison.cols <- outer( north.cells, north.cells, function (x,y) { ifelse( xor(x,y), "red", ifelse(x, "blue", "purple") ) } )

    # actually run ms
    run_ms( dem, nsamp=sample.config, trees=TRUE, outdir=outdir, nreps=ntrees )

    tree.output <- trees_from_ms( file.path(outdir,"msoutput.txt") )
    tree.dists <- tree_dists( tree.output, sample.config=sample.config )
    mean.dist <- tree.dists[[1]]
    for (k in seq_along(tree.dists)[-1]) { mean.dist <- mean.dist+tree.dists[[k]] }
    mean.dist <- mean.dist/length(tree.dists)

    # make things necessary to make nice plots of the results
    sample.ids <- paste0("sample-",1:sum(sample.config[,"n"]))
    ut <- upper.tri(mean.dist,diag=FALSE)
    dist.table <- data.frame( id1=sample.ids[row(mean.dist)[ut]], id2=sample.ids[col(mean.dist)[ut]], pi=mean.dist[ut] )
    geog.table <- data.frame( id1=sample.ids[row(mean.dist)[ut]], id2=sample.ids[col(mean.dist)[ut]], distance=geog.dist[l2s,l2s][ut] )
    sample.spatialpoints <- sample.points[l2s]
    row.names(sample.spatialpoints) <- sample.ids
    fake.pcs <- data.frame( id=sample.ids, pc1=ifelse(north.cells[l2s],+1.0,-1.0) )
    # write out things necessary to make nice plots of the results
    write.csv( dist.table, file=file.path(outdir,"mean-divergence.csv"), row.names=FALSE )
    write.csv( geog.table, file=file.path(outdir,"geog-distance.csv"), row.names=FALSE )
    save( sample.spatialpoints, file=file.path(outdir,"geog_coords.RData"))
    write.csv( fake.pcs, file=file.path(outdir,"pcs.csv"), row.names=FALSE )

    # actual pcs
    #  of the covariance matrix (up to scaline; see McVean)
    covmat <- (rowMeans(mean.dist) + colMeans(mean.dist) - mean(mean.dist) - mean.dist)
    pmat <- diag(nrow(covmat)) - 1/nrow(covmat)
    pcs <- eigen( (pmat %*% covmat %*% t(pmat)) )$vectors[,1:5]
    pc.pal <- rainbow(n=32, start=4/6, end=0)
    pc.cols <- apply(pcs, 2, function (x) { pc.pal[cut(x,length(pc.pal))] } )

    ## look at some plots

    pdf(file=file.path(outdir,"trees-and-things.pdf"),width=10,height=5,pointsize=10)

    plot( as.vector(geog.dist[l2s,l2s])/1e3, as.vector(mean.dist), pch=20, cex=0.5, 
         xlab="geog dist (km)", ylab="mean TMRCA (y)",
         col=adjustcolor(comparison.cols[l2s,l2s],0.5) )

    for (k in 1:ncol(pcs)) {
        plot( crop(habitat,sample.xy), main=sprintf("PC %d",k) )
        points( sample.xy, pch=20, cex=2, col=pc.cols[,k] )
        # if (interactive() && is.null(locator(1))) { break }
    }

    layout(t(1:2))
    for (tree in tree.output[1:plot.ntrees]) {
        plot( crop(habitat,sample.xy) )
        points( sample.xy, pch=20, cex=2, col=pc.cols[,1] )
        # plot_sample_config( dem, sample.config, add=TRUE, xy=pop.xy, col=sample.cols )
        pp <- plot.phylo( tree, show.tip.label=FALSE )  # tip.color=pc.cols[,k] )  # tip.color=sample.cols[l2s], cex=0.2 )
        axisPhylo(1)
        ab <- abline_phylo(v=dem@t[length(dem@t)], lty=2, col='grey')
        points( rep(1.05*ab[1],nrow(pcs)), 1:nrow(pcs), pch=20, col=pc.cols[,1][order(tfn(tree))] )
        # if (interactive() && is.null(locator(1))) { break }
    }

    dev.off()
}
