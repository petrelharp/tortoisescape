#!/usr/bin/Rscript

usage <- "
     Fit parameters to a set of landscape layers,
         Rscript analytic-parameter-inference.R (input json configuration file) (output json file)
       e.g.
         Rscript analytic-parameter-inference.R test_six_layers/noise-0.01/256x/config.json test_six_layers/noise-0.01/256x/inference.json

"
if ( !interactive() && length(commandArgs(TRUE))<2 ) { stop(usage) }

source("input-output-fns.R")
opts <- read.config("config"="json","outfile"="text")
for (x in config$setup_files) {
    loadfile <- file.path( dirname(opts[['config']]), x )
    cat("Loading ", x, " .\n")
    load( loadfile )
}

# this will store the output
output.config <- config

# match these hitting times
obs.ht <- read.sub.hts( config$observed_ht_file, locs )

source("resistance-fns.R")
require(parallel)
numcores<-getcores()

for (this.step in 1:config$maxstep) {
    init.params <- paramvec(output.config)

    # Set up inference
    G@x <- update.G( init.params )
    dG <- rowSums(G)

    # interpolate hitting times based on current parameters
    hts <- interp.hitting( neighborhoods, G-diag(rowSums(G)), obs.ht, obs.locs=locs, alpha=config$alpha, numcores=numcores )

    # infer parameters based on full hitting times
    zeros <- unlist(neighborhoods) + rep((seq_along(neighborhoods)-1)*nrow(hts),sapply(neighborhoods,length))
    scaling <- 1 # sqrt(nrow(G) * length(locs))
    sc.one <- 1/scaling
    hts <- hts/scaling
    hts[zeros] <- 0

    LdL <- params.logistic.setup()
    L <- LdL$L
    dL <- LdL$dL

    parscale <- unlist( config$paramscale )
    Lval <- L( init.params )

    param.optim <- optim( par=init.params, fn=L, gr=dL, control=list(parscale=parscale,fnscale=max(1,abs(Lval)/10)), method="BFGS" )
    if (param.optim$convegence != 0) { warning("optim() failed to converge.  Continuing anyhow.") }

    paramvec(output.config) <- param.optim$par

    if (all( abs(param.optim$par-paramvec(output.config)) < unlist(config$paramtol) )) { 
        output.config$optim.results <- param.optim
        output.config$optim.nsteps <- this.step
        break
    }
}

if (this.step==config$maxstep) { warning("Maximum number of steps reached.") }

write.json.config( output.config, outfile )

