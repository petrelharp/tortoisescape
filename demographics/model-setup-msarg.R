source("msarg/msarg.R")
library(raster)
library(landsim)
library(sp)
library(rgeos)

# load tortoise habitat
habitat <- raster("../visualization/nussear_masked.grd")
sub.hab <- aggregate( habitat, fact=32 )

if (FALSE) {
    # delineate the 'northern' area
    plot(habitat)
    norths <- drawPoly()
    proj4string(norths) <- CRS(proj4string(habitat))
    save(norths,file="../visualization/north-polygon.RData")
} else {
    load("../visualization/north-polygon.RData")
}

###
# PARAMETERS
# we'll construct a migration matrix across scaling.fac generations to get migration between coarse cells
tort.sigma <- 0.5e3 # m/gen
scaling.fac <- floor(mean(res(sub.hab)/tort.sigma)^2)/5  # generations

# set up things
pop <- population( sub.hab, genotypes='a', accessible=rep(TRUE,length(values(sub.hab))), 
                  habitable=( !is.na(values(sub.hab)) & values(sub.hab)>0 ) )
pop.xy <- xyFromCell(sub.hab,cell=which(pop$accessible))
pop.points <- SpatialPoints(pop.xy, proj4string=CRS(proj4string(habitat)))
# only nearest-neighbor migration (king's neighborhood)
migr <- migration( kern="gaussian", sigma=sqrt(scaling.fac)*tort.sigma, 
                  radius=res(sub.hab), normalize=1 )

# create a population object for msarg
npops <- sum(pop$accessible)
# density: one tortoise per 20 ha in perfect habitat
tort.density <- ( 0.05 * prod(res(pop$habitat)) / 1e4 ) * values(pop$habitat)[pop$accessible]
tort.density[is.na(tort.density)] <- 0
migr.mat <- migration_matrix( pop, migration=migr ) / scaling.fac
diag(migr.mat) <- 0
hab.grid <- new("gridArray",
                npop = as.integer(dim(pop$habitat)[c(2,1,3)]),
                N = tort.density,
                G = rep( 0, npops ),
                M = migr.mat
            )

# choose where to sample things at: 2 in each location
nsamps <- 1  # per location
sample.cells <- which( values(sub.hab) > 0.5 )
sample.config <- arrayInd( sample.cells, .dim=dim(hab.grid) )
colnames(sample.config) <- c("row","col","layer")
sample.config <- sort_sample_config( cbind(sample.config,n=nsamps) )

plot(habitat)
plot(hab.grid,xy=pop.xy,add=TRUE)

###
# create a barrier
barrier.end.time <- 20e3/25    # 20 Kya
barrier.begin.time <- 500e3/25  # 200 Kya

barrier.demog <- demography( hab.grid )

north.cells <- gContains( norths, pop.points, byid=TRUE )
barrier.mat <- migr.mat
barrier.mat@x[ xor(north.cells[rowinds(barrier.mat)],north.cells[colinds(barrier.mat)]) ] <- 0.0

barrier.demog <- add_to_demography( barrier.demog, dt=barrier.end.time )
barrier.demog[[2]]@M <- barrier.mat

###
# end the barrier

barrier.demog <- add_to_demography( barrier.demog, tnew=barrier.begin.time, pop=barrier.demog[[1]] )

if (!file.exists("barrier-demography-msarg.RData")) {
    save(barrier.demog, sample.config, pop.xy, file="barrier-demography-msarg.RData")
} else {
    warning("File 'barrier-demography-msarg.RData' already exists.")
}

plot(barrier.demog,xy=pop.xy)


###
# create the refugia
expansion.end.time <- 1e4/25
refugia.end.time <- 20e4/25
refugia.begin.time <- 50e4/25
expansion.speed <- 3e3/mean(res(pop$habitat))/(refugia.end.time-expansion.end.time)
expansion.width <- 100e3/mean(res(pop$habitat))

refugia.centers <- SpatialPoints(coords=cbind( x=c(727229.8,626639.7),
                                               y=c(3803747,3957436)),
                        proj4string=CRS(proj4string(habitat)))
refugia.radius <- 5e4
refugia <- gBuffer( refugia.centers, width=refugia.radius, byid=TRUE )
refugia.mask <- gContains( gUnaryUnion(refugia), pop.points, byid=TRUE )

refugia.demog <- demography( hab.grid )
refugia.demog <- add_to_demography( refugia.demog, tnew=refugia.end.time,
                                   fn=modify_grid_layer, layer=1, dN=refugia.mask )
refugia.demog <- add_to_demography( refugia.demog, tnew=refugia.begin.time,
                                   pop=refugia.demog[[1]] )
refugia.demog <- logistic_interpolation( refugia.demog, 
                                         t.end=expansion.end.time,
                                         t.begin=refugia.end.time, 
                                         nsteps=30, 
                                         speed=expansion.speed, 
                                         width=expansion.width )
if (!check_demography(refugia.demog)) { stop("Final state not communicating.") }

matplot(t(sapply(as(refugia.demog,"list"),slot,"N")),type='l')

plot(refugia.demog,xy=pop.xy)

# check for jumps that are too big that things get stuck
for (k in 2:length(refugia.demog)) {
    dothese <- xor( refugia.demog[[k]]@N==0, refugia.demog[[k-1]]@N==0 )
    if (any(dothese)) {
        plot( pop.xy[ dothese, ], pch=20, cex=2, 
             col=ifelse(refugia.demog[[k]]@N==0,"black","red")[dothese],
            main=sprintf("t=%f",refugia.demog@t[k-1]) )
        if (interactive() && is.null(locator(1))) { break }
    }
}

if (!file.exists("refugia-demography-msarg.RData")) {
    save(refugia.demog, sample.config, pop.xy, file="refugia-demography-msarg.RData")
} else {
    warning("File 'refugia-demography-msarg.RData' already exists.")
}

