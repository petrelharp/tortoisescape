#' The functions here should take the demographic parameters (see README.md),
#' set up a demographic model, simulate data,
#' and compare the results to the observed data.


#' Set up the demographic model
#'
#' @param pop A population object.
#' @param pop.density Population density, in indivs/m^2
#' @param sigma Migration rate, in meters/generation.
#' @param refugia.coords Centers of refugia, in UTM coordinates (meters).
#' @param refugia.radii Common radius of refugia, in meters.
#' @param refugia.time Time refugia began, in generations.
#' @param expansion.time Time refugia ended and expansion began, in generations.
#' @param expansion.speed Speed of expansion, in meters per generation.
#' @param expansion.width Width of the expansion, in meters.
model_setup <- function (
                         pop,
                         pop.density,
                         sigma,
                         refugia.coords,
                         refugia.radii,
                         refugia.time,
                         contraction.speed,
                         expansion.time,
                         expansion.speed,
                         expansion.width
                     ) {

    ######
    # landsim pre-setup
    pop.xy <- xyFromCell(pop$habitat,cell=which(pop$accessible))
    pop.points <- SpatialPoints(pop.xy, proj4string=CRS(proj4string(pop$habitat)))

    # msarg pre-setup
    density.vec <- ( pop.density * prod(res(pop$habitat)) ) * values(pop$habitat)[pop$accessible]
    density.vec[is.na(density.vec)] <- 0

    # create the base migration matrix: WILL BE MODIFIED
    migr <- migration( kern="gaussian", sigma=sigma, 
                      radius=3*sigma, normalize=1, discretize=TRUE,
                      disc.fact=3*mean(res(pop$habitat))/sigma )
    migr.mat <- migration_matrix( pop, migration=migr )
    diag(migr.mat) <- 0


    # create a population object for msarg
    pop.N <- ( pop.density * prod(res(pop$habitat)) ) * values(pop$habitat)[pop$accessible]
    pop.N[is.na(pop.N)] <- 0
    hab.grid <- new("gridArray",
                    npop = as.integer(dim(pop$habitat)[c(2,1,3)]),
                    N = pop.N,
                    G = rep( 0, sum(pop$accessible) ),
                    M = migr.mat
                )

    # create the refugia  WILL BE MODIFIED
    refugia.centers <- SpatialPoints(coords=refugia.coords,
                            proj4string=CRS(proj4string(pop$habitat)))
    refugia <- gBuffer( refugia.centers, width=refugia.radii, byid=TRUE )
    refugia.mask <- gContains( gUnaryUnion(refugia), pop.points, byid=TRUE )
    refugia.demog <- demography( hab.grid )
    # endpoint of refugia
    refugia.demog <- add_to_demography( refugia.demog, tnew=expansion.time,
                                       fn=modify_grid_layer, layer=1, dN=refugia.mask )
    # beginning of refugia: before was the same as now (abruptly)
    refugia.demog <- add_to_demography( refugia.demog, tnew=refugia.time,
                                       pop=refugia.demog[[1]] )
    # add re-expansion to modern day (t.end=0):
    #  note due to speed may effectively end much earlier
    refugia.demog <- logistic_interpolation( refugia.demog, 
                                             t.end=0,
                                             t.begin=expansion.time, 
                                             nsteps=30, 
                                             speed=expansion.speed/mean(res(pop$habitat)), 
                                             width=expansion.width/mean(res(pop$habitat)) )
    # check this works
    if (!check_demography(refugia.demog)) { stop("Final state not communicating.") }

    return( refugia.demog )
}

#' Run ms and compute mean pairwise distances
#' @param dem A demography object, as produced by model_setup().
#' @param sample.config A sample configuration, as desired by run_ms().
#' @param outdir The output directory.
#' @param ntrees The number of trees to simulate.
sim_data <- function (
                      dem,
                      sample.config,
                      outdir,
                      ntrees
                      ) {
    # actually run ms
    run_ms( dem, nsamp=sample.config, trees=TRUE, outdir=outdir, nreps=ntrees, N.eps=0.1 )

    # parse ms output
    tree.output <- trees_from_ms( file.path(outdir,"msoutput.txt") )
    tree.dists <- tree_dists( tree.output, sample.config=sample.config )
    mean.dist <- tree.dists[[1]]
    for (k in seq_along(tree.dists)[-1]) { mean.dist <- mean.dist+tree.dists[[k]] }
    mean.dist <- mean.dist/length(tree.dists)

    return( mean.dist )
}


#' Makes a plot like `trees-and-things.pdf`.
plot_results <- function ( 
                          ) {
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
