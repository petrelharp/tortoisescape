#!/usr/bin/Rscript

usage <- "
Infer parameters, using the (slow) direct method on a single resolution.
Usage:
    Rscript direct-inference.R (config file) (output file) (max steps) [result from previous run]
e.g.
    Rscript sim-hitting-times.R real-data/six-layer-256x.json real-data/six-layer-256x-trust_1.RData 100

"

argvec <- if (!interactive()) { commandArgs(TRUE) } else { scan(what='char') }

if (length(argvec)<3) { stop(usage) }
config.file <- argvec[1]
output.file <- argvec[2]
maxit <- as.numeric( argvec[3] )
prev.file <- if (length(argvec)>3) { argvec[4] } else { NULL }

source("resistance-fns.R")
require(parallel)
require(trust)

config <- read.json.config(config.file)
layer.names <- config$layer_names
for (x in file.path(dirname(config.file),config$setup_files)) {
    load( x )
}

# do relative to reference individuals
ref.inds <- which.nonoverlapping(neighborhoods)

if (is.null(prev.file)) {
    init.params <- paramvec(config)
} else {
    load(prev.file)
    init.params <- trust.optim$argument
}

G@x <- update.G(init.params[-1])

# the setup
ds <- direct.setup( obs.locs=locs, obs.hts=pimat[,ref.inds], 
        neighborhoods=neighborhoods[ref.inds], 
        G=G, update.G=update.G, layers=layers, 
        transfn=transfn, valfn=valfn, ndelta=ndelta, ngamma=ngamma
    )

# check this runs
system.time( init.value <- ds(init.params) )
init.value$value  

trust.optim <- trust( objfun=ds, parinit=init.params, rinit=0.25, rmax=5, iterlim=maxit, blather=TRUE )

trust.optim$ref.inds <- ref.inds
trust.optim$config.file <- config.file
trust.optim$invocation <- paste(commandArgs())
trust.optim$prev.file <- prev.file

save( trust.optim, file=output.file )
