library(raster)
library(sp)
library(rgeos)
library(ape)
library(jsonlite)
library(msarg)
library(landsim)

# load tortoise habitat
hab.fact <- 32
full.habitat <- raster("../visualization/nussear_masked.grd")
habitat <- aggregate( full.habitat, fact=hab.fact )

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

good.params <- read.csv("good-results-lt-1e6.csv",header=TRUE,stringsAsFactors=FALSE)
init.params.vec <- as.numeric( good.params[ sample.int(nrow(good.params),1), sapply(good.params,is.numeric)] )
names(init.params.vec) <- names(good.params)[sapply(good.params,is.numeric)]
init.params <- list(
        pop.density = init.params.vec["pop.density"],       # density in indivs/m^2 
        sigma = init.params.vec["sigma"],             # m/gen
        refugia.coords = cbind( x=init.params.vec[c("refugia.coords.1","refugia.coords.2")],
                                y=init.params.vec[c("refugia.coords.3","refugia.coords.4")]),   # in meters (UTM)
        refugia.radii = init.params.vec["refugia.radii"],               # in meters
        refugia.time = init.params.vec["refugia.time"],      # length of time refugia existed for, in generations
        expansion.time = init.params.vec["expansion.time"],    # time refugia ended, in generations
        expansion.speed = init.params.vec["expansion.speed"],  # in m/gen (not really??)
        expansion.width = init.params.vec["expansion.width"]   # in meters, roughly
)
# bounds on reasonable places (?) for refugia to be centered
refugia.bounds <- list( x=c( 350000, 780000),
                        y=c(3640000,4200000) )

### Other (unchanging) parameters
ntrees <- 500    # Number of trees to simulate.

run.id <- floor(1e6*runif(1))
basedir <- sprintf("run_%06d",run.id)
dir.create( basedir )
cat("Writing out to", basedir, "----\n")

cat( toJSON(c(init.params,list(ntrees=ntrees,hab.fact=hab.fact)), pretty=TRUE), file=file.path(basedir,"params.json") )
write.csv( dist.df, file=file.path(basedir,"dist-df.csv"), row.names=FALSE )

# explore randomly nearby with these
n.iter <- 100

which.max.score <- 1
max.score <- Inf
# explore parameter space
# iter.num <- 1
for (iter.num in seq_len(n.iter)) {

    if (iter.num == 1) {
        params <- init.params
    } else { 
        params <- lapply( init.params, function (x) {
                    x * (0.5+runif(length(x)))
            } )
        params$refugia.coords <- sampleRandom(habitat,size=2,xy=TRUE)[,1:2]
    }

    score <- run_sim( params, iter.num, ntrees )
    if (score < max.score) { max.score <- score; init.params <- params; which.max.score <- iter.num }
}

# to translate parameter vector to params
param.names <- names(init.params)
param.inds <- lapply( seq_along(init.params), function (k) {
            sum( unlist(sapply( init.params[ seq_len(k-1) ], length )) ) + seq_len(length(init.params[[k]]))
        } )
optim_fun <- function (x) {
    # put parameters in correct spots
    params <- init.params
    for (k in seq_along(params)) {
        params[[k]][] <- x[param.inds[[k]]]
    }
    out <- run_sim( params, iter.num=floor( 1e6*runif(1) ), ntrees=100 )
    cat("score:", out,"\n")
    return(out)
}

# otherwise do optimization
optim( unlist(init.params), optim_fun, 
     lower=c(0,0,
             refugia.bounds[['x']][1],refugia.bounds[['x']][1],
             refugia.bounds[['y']][1],refugia.bounds[['y']][1],
             0,0,0,0,0), 
     upper=c(Inf,Inf,
             refugia.bounds[['x']][2],refugia.bounds[['x']][2],
             refugia.bounds[['y']][2],refugia.bounds[['y']][2],
             Inf,Inf,Inf,Inf,Inf), 
     method="L-BFGS-B",
     control=list( 
                 fnscale=1e4,
                 parscale=c(1e-06, 1e+05, 1e+06, 1e+06, 1e+07, 1e+07, 1e+05, 1e+05, 1e+05, 1e-00, 1e+05),
                 maxit=25 )
            )
