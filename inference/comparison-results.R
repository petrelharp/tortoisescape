#!usr/bin/Rscript

usage <- "Use a fitted model to make predictions in a way comparable to others.
Usage:
    Rscript comparison-results.R (name of results .RData file) [name of output .RData file]
and the output defaults to the input file, but with -comparison.RData appended.
"

argvec <- if (interactive()) { scan(what='char') } else { commandArgs(TRUE) }
if (length(argvec)<1) { stop(usage) }

source("resistance-fns.R")
require(parallel)
numcores<-getcores()
require(raster)
require(rgdal)

infile <- argvec[1]
load(infile)
outfile <- if (length(argvec)>1) { argvec[2] } else { gsub("[.]R[dD]ata$", "-comparison.RData", infile) }

# setup for using *everyone*, on habitat only
config.file <- "summaries/habitat-only/config.json"
config <- read.json.config(config.file)
for (x in config$setup_files) { load(file.path(dirname(config.file),x)) }
config$reference_inds <- row.names( sample.locs )

# and for this model
local.config <- read.json.config( trust.optim$config.file )
layer.names <- local.config$layer_names
paramvec(local.config) <- trust.optim$argument

# get locally-specified layers in globally-specified way
use.files <- layer.files[match(layer.names,layer.file.names)]
names(use.files) <- layer.names

## find if other layers have additional NA values, and block these later
new.nalayer <- nalayer
for (lf in use.files) {
    other <- raster(lf) 
    new.nalayer <- mask( new.nalayer, other )
}
block.these <- match( which( is.na(values(new.nalayer)) & ! is.na(values(nalayer)) ), nonmissing )

# load up the layers
if (length(layer.names)>0) { 
    layers <- sapply( layer.names, function (ln) {
                scale( values( raster( paste(full.layer.prefix,ln,sep='') ) )[nonmissing] )
            } )
} else {
    layers <- matrix(0,nrow=nrow(G),ncol=0)
}
stopifnot(nrow(layers)==nrow(G))
# ADD the constant layer
layers <- cbind( 1, layers )
layer.names <- c( "constant", layer.names )

G@x <- update.G(paramvec(local.config)[-1])

hts <- hitting.analytic( neighborhoods, G, numcores=numcores, blocked=block.these )

trust.optim.results <- trust.optim[ "converged" ]
    
save( infile, config, local.config, trust.optim.results, hts, pimat, file=outfile )
