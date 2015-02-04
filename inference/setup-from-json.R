#!/usr/bin/Rscript

usage <- "Set up for inference, from a json config file.
Usage:
    Rscript setup-from-json.R (name of json file) (output RData file)
e.g.
    Rscript setup-from-json.R real-data/nussear-256x/config.json real-data/nussear-256x/setup.RData
"

argvec <- if (interactive()) { scan(what='char') } else { commandArgs(TRUE) }
if (length(argvec)<1) { stop(usage) }

source("resistance-fns.R")
# require(parallel)
numcores<-getcores()
require(raster)
require(rgdal)

config.file <- argvec[1]
config <- read.json.config(config.file)
outfile <- if (length(argvec)>1) { argvec[2] } else { file.path(dirname(config.file),config$setup_files) }

# save everything after this point
source.ls <- ls()

layer.names <- config$layer_names
full.layer.prefix <- file.path(dirname(config.file),config$layer_prefix)
dir.layer.prefix <- if (grepl("/$",config$layer_prefix)) { full.layer.prefix } else { dirname(full.layer.prefix) }
base.layer.prefix <- if (grepl("/$",config$layer_prefix)) { "" } else { basename(config$layer_prefix) }

### from make-overlap-na-layer.R

# setup NA overlap layer

layer.files <- gsub("//","/",list.files(full.layer.prefix,pattern="gri$",full.names=TRUE),fixed=TRUE)
layer.file.names <- gsub(paste(".*",config$layer_prefix,sep=''),"",gsub("[.]((gri)|(tif)|(asc))$","",layer.files))
use.files <- layer.files[match(layer.names,layer.file.names)]
names(use.files) <- layer.names

nalayer <- raster( grep( paste(config$mask_layer,"[.]((gri)|(tif)|(asc))$",sep=''), layer.files, value=TRUE )[1] )

for (lf in use.files) {
    other <- raster(lf) 
    nalayer <- mask( nalayer, other )
}


### from setup-real-G.R

nonmissing <- which(!is.na(values(nalayer)))

G <- make.G( layer=nalayer, nonmissing=nonmissing )
Gjj <- p.to.j(G@p)

if (length(layer.names)>0) { 
    layers <- sapply( layer.names, function (ln) {
                values( raster( paste(full.layer.prefix,ln,sep='') ) )[nonmissing]
            } )
    layer.center <- apply(layers,2,mean,na.rm=TRUE)
    layer.scale <- apply(layers,2,sd,na.rm=TRUE)
    layers <- sweep(sweep(layers,2,layer.center,"-"),2,layer.scale,"/")
} else {
    layers <- matrix(0,nrow=nrow(G),ncol=0)
    layer.scale <- layer.center <- numeric(0)
}
stopifnot(nrow(layers)==nrow(G))

# ADD the constant layer
layers <- cbind( 1, layers )
layer.center <- c(0,layer.center)
layer.scale <- c(1,layer.scale)
layer.names <- names(layer.center) <- names(layer.scale) <- c( "constant", layer.names )

# transfn <- exp
transfn <- function (x) { 1/(1+exp(-x)) }
# valfn <- function (gamma) { ( rowSums( layers * gamma[col(layers)], na.rm=TRUE ) ) }
# this is faster if we don't have to worry about NAs (we shouldn't?)
valfn <- function (gamma) { ans <- layers[,1]*gamma[1]; for (k in (1:NCOL(layers))[-1]) { ans <- ans+layers[,k]*gamma[k] }; return(ans) }

ndelta <- ngamma <- length(layer.names)
update.G <- function(params) {
    stopifnot(length(params)==1+2*ncol(layers))
    beta <- params[1]
    gamma <- params[1+(1:ngamma)]
    delta <- params[1+ngamma+(1:ndelta)]
    return( exp(beta) * transfn(valfn(gamma))[G@i+1L] * transfn( valfn(delta)[G@i+1L] + valfn(delta)[Gjj] ) )
}

G@x <- update.G(paramvec(config)[-1])


### from setup-tort-locs.R

sample.loc.obj <- load(file.path(dirname(config.file),config$sample_locs))
assign("sample.locs",get(sample.loc.obj))
sample.locs <- spTransform( sample.locs, CRSobj=CRS(proj4string(nalayer)))
rm(sample.loc.obj)
sample.ids <- row.names(sample.locs)

orig.locs <- cellFromXY( nalayer, sample.locs )
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
