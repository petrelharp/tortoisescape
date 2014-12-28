#!/usr/bin/Rscript

usage <- "
     Do MCMC over the parameters, matching hitting times to observed values.
         Rscript landscape-mcmc.R (input json configuration file) (output json config file) (output MCMC trace file)
       e.g.
         Rscript landscape-mcmc.R test_six_layers/noise-0.00/256x/true-config.json test_six_layers/noise-0.00/256x/mcmc-results.json test_six_layers/noise-0.00/256x/mcmc-trace.tsv

"
if ( !interactive() && length(commandArgs(TRUE))<2 ) { stop(usage) }

source("input-output-fns.R")
opts <- read.config("config"="json","outfile"="text","outtrace"="text")
this.seed <- getset.seed()
# this will store output things
output.config <- config
output.config$seed <- this.seed

source("resistance-fns.R")
numcores <- getcores()

# file names in config are relative to whereever it lives
config.dir <- dirname(opts[['config']])

# load input setup files
for (x in config$setup_files) {
    loadfile <- file.path( config.dir, x )
    cat("Loading ", x, " .\n")
    load( loadfile )
}

# match these hitting times
obs.ht <- read.sub.hts( file.path(config.dir,config$observed_ht_file), locs )
parscale <- unlist( config$paramscale )

M <- function (params) {
    G@x <- update.G(params)
    hts <- hitting.analytic(neighborhoods, G-diag(rowSums(G)), numcores=numcores )
    ans <- (-1) * ( sum( ( hts[locs,] - obs.ht )^2 ) )
    return( if (!is.finite(ans)) { Inf } else { ans } )
}

require(mcmc)

mcmc.run <- metrop( M, initial=paramvec(config), nbatch=100, blen=1, scale=parscale )
