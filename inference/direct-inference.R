#!/usr/bin/Rscript

usage <- "
Infer parameters, using the (slow) direct method on a single resolution.
Usage:
    Rscript direct-inference.R (config file) (output file) (max steps) [result from previous run]
e.g.
    Rscript sim-hitting-times.R real-data/six-layer-256x.json real-data/six-layer-256x-trust_1.RData 100

"

argvec <- if (!interactive()) { commandArgs(TRUE) } else { scan(what='char') }
cat("argvec:\n")
cat(paste(argvec,collapse=" "),"\n")

if (length(argvec)<3) { stop(usage) }
config.file <- argvec[1]
output.file <- argvec[2]
maxit <- as.numeric( argvec[3] )
prev.file <- if (length(argvec)>3) { argvec[4] } else { NULL }

source("resistance-fns.R")
require(parallel)
require(trust)

cat("Reading config: ", config.file, "\n")
config <- read.json.config(config.file)
layer.names <- config$layer_names
for (x in file.path(dirname(config.file),config$setup_files)) {
    cat("Loading: ", x, "\n")
    load( x )
}

# do relative to reference individuals
ref.inds <- if (is.null(config$reference_inds)) {
        which.nonoverlapping(neighborhoods)
    } else {
        match( config$reference_inds, rownames(pimat) )
    }
if (any(is.na(ref.inds))) {
    warning(paste("Removed", paste(setdiff(config$reference_inds,rownames(pimat)),collapse=", "), "from reference individuals."))
}
ref.inds <- ref.inds[ !is.na(ref.inds) ]

if (is.null(prev.file)) {
    init.params <- paramvec(config)
} else {
    cat("Restarting from ", prev.file, "\n")
    load(prev.file)
    init.params <- trust.optim$argument
}
param.scale <- paramvec(config,"param_scale")

G@x <- update.G(init.params[-1])

# the setup
#   omit self comparisons from the fitting procedure
fit.pimat <- pimat[ref.inds,ref.inds]  
# note:  changing to symmetric=FALSE can:
#   remove subsetting to rows if also set symmetric=FALSE below
#   and change obs.locs to locs
obs.locs <- locs[ref.inds]
diag(fit.pimat) <- NA
ds <- direct.setup( obs.locs=obs.locs, obs.hts=fit.pimat,
        neighborhoods=neighborhoods[ref.inds], 
        G=G, update.G=update.G, layers=layers, 
        transfn=transfn, valfn=valfn, ndelta=ndelta, ngamma=ngamma,
        symmetric=TRUE
    )

# check this runs
system.time( init.value <- ds(init.params) )
init.value$value  

cat("Running trust().\n")

# 'parscale' according to 'trust' should be 1/(typical step size)
trust.optim <- trust( objfun=ds, parinit=init.params, parscale=1/param.scale,
        rinit=0.25, rmax=5, iterlim=maxit, blather=TRUE )

cat("Done with trust().\n")

trust.optim$param.scale <- param.scale
trust.optim$ref.inds <- ref.inds
trust.optim$config.file <- config.file
trust.optim$invocation <- paste(commandArgs())
trust.optim$prev.file <- prev.file
trust.optim$pimat <- fit.pimat

cat("Saving to ", output.file,"\n")
save( trust.optim, file=output.file )
