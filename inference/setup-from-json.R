#!/usr/bin/Rscript

usage <- "Set up for inference, from a json config file.
Usage:
    Rscript setup-from-json.R (name of json file) (output RData file)
e.g.
    Rscript setup-from-json.R real-data/nussear-256x/config.json real-data/nussear-256x/setup.RData
"

argvec <- if (interactive()) { scan(what='char') } else { commandArgs(TRUE) }
if (length(argvec)<2) { stop(usage) }

source("resistance-fns.R")
# require(parallel)
numcores<-getcores()
require(raster)
rasterOptions(tmpdir=".")

config.file <- argvec[1]
outfile <- argvec[2]
config <- read.json.config(config.file)

source.ls <- ls()

layer.names <- config$layer_names
full.layer.prefix <- file.path(dirname(config.file),config$layer_prefix)

### from make-overlap-na-layer.R

# setup NA overlap layer

layer.files <- list.files(dirname(full.layer.prefix),paste(basename(config$layer_prefix),".*gri",sep=''),full.names=TRUE)
layer.file.names <- gsub(paste(".*",config$layer_prefix,sep=''),"",gsub(".gri","",layer.files,fixed=TRUE))
use.files <- layer.files[match(layer.names,layer.file.names)]
names(use.files) <- layer.names

# refname <- grep("dem_30",layer.file.names)[1]
# nalayer <- is.na( raster(layer.files[refname]) ) # dem
na.layer.base <- if ( is.null(config$mask.layer) ) { "dem_30" } else { config$mask.layer }
refname <- grep( na.layer.base, layer.file.names )[1]
nalayer <- is.na( raster(config$mask.layer) )

for (lf in use.files) {
    other <- raster(lf) 
    nalayer <- nalayer | is.na(other)
}


### from setup-real-G.R

nonmissing <- which(values(nalayer)==0)

G <- make.G( layer=nalayer, nonmissing=nonmissing )
Gjj <- p.to.j(G@p)

layers <- sapply( layer.names, function (ln) {
            scale( values( raster( paste(full.layer.prefix,ln,sep='') ) )[nonmissing] )
        } )
stopifnot(nrow(layers)==nrow(G))

# ADD the constant layer
layers <- cbind( 1, layers )
layer.names <- c( "constant", layer.names )

# transfn <- exp
transfn <- function (x) { 1/(1+exp(-x)) }
# valfn <- function (gamma) { ( rowSums( layers * gamma[col(layers)], na.rm=TRUE ) ) }
# this is faster if we don't have to worry about NAs (we shouldn't?)
valfn <- function (gamma) { ans <- layers[,1]*gamma[1]; for (k in (1:NCOL(layers))[-1]) { ans <- ans+layers[,k]*gamma[k] }; return(ans) }

ndelta <- ngamma <- length(layer.names)
update.G <- function(params) {
    beta <- params[1]
    gamma <- params[1+(1:ngamma)]
    delta <- params[1+ngamma+(1:ndelta)]
    return( exp(beta) * transfn(valfn(gamma))[G@i+1L] * transfn( valfn(delta)[G@i+1L] + valfn(delta)[Gjj] ) )
}

G@x <- update.G(paramvec(config))


### from setup-tort-locs.R

sample.loc.obj <- load(file.path(dirname(config.file),config$sample_locs))
assign("sample.locs",get(sample.loc.obj))
rm(sample.loc.obj)
sample.ids <- row.names(sample.locs)

orig.locs <- cellFromXY( nalayer, tort.coords.rasterGCS )
locs <- match(orig.locs,nonmissing)

ndist <- 15000  # 15 km

neighborhoods <- get.neighborhoods( ndist=ndist, locations=sample.locs, nonmissing=nonmissing, layer=nalayer, numcores=numcores, na.rm=FALSE )

### from setup-inference.R

# load in 
pimat.df <- read.csv( file.path(dirname(config.file),config$divergence_file) )
pimat <- .018 * 1e8 * df.to.sym(pimat.df,sample.ids)
stopifnot(sum(is.na(pimat))==0)

# remove NA individuals
na.indiv <- is.na( locs )
locs <- locs[!na.indiv]
sample.ids <- sample.ids[!na.indiv]
sample.locs <- sample.locs[!na.indiv]
orig.locs <- orig.locs[!na.indiv]
neighborhoods <- neighborhoods[!na.indiv]
pimat <- pimat[!na.indiv,!na.indiv]

save(list=setdiff(ls(),source.ls), file=outfile)
