#!/usr/bin/Rscript

usage <- "
     Fit parameters to a set of landscape layers,
         Rscript analytic-parameter-inference.R (input json configuration file) (output json file)
       e.g.
         Rscript analytic-parameter-inference.R test_six_layers/noise-0.00/256x/true-config.json test_six_layers/noise-0.00/256x/true-inference.json

"
if ( !interactive() && length(commandArgs(TRUE))<2 ) { stop(usage) }

source("input-output-fns.R")
opts <- read.config("config"="json","outfile"="text")
this.seed <- getset.seed()
# this will store the output
output.config <- config
output.config$seed <- this.seed

source("resistance-fns.R")

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

# setup
zeros <- unlist(neighborhoods) + rep((seq_along(neighborhoods)-1)*nrow(G),sapply(neighborhoods,length))
parscale <- unlist( config$paramscale )

require(parallel)
numcores<-getcores()

for (this.step in 1:config$maxstep) {

    init.params <- paramvec(output.config)
    # Set up inference
    G@x <- update.G( init.params )
    dG <- rowSums(G)

    # interpolate hitting times based on current parameters
    hts <- interp.hitting( neighborhoods, G-diag(dG), obs.ht, obs.locs=locs, alpha=config$alpha, numcores=numcores )
    if (!is.numeric(hts)) { print.and.dump() }

    # infer parameters based on full hitting times
    sc.one <- 1
    hts[zeros] <- 0

    LdL <- params.logistic.setup(init.params,G,update.G,hts,zeros,sc.one,layers,transfn,valfn,ndelta,ngamma)
    L <- LdL$L
    dL <- LdL$dL

    Lval <- L( init.params )

    param.optim <- optim( par=init.params, fn=L, gr=dL, control=list(parscale=parscale,fnscale=max(1,abs(Lval)/10),maxit=config$maxit), method="BFGS" )
    converged <- (param.optim$convergence == 0)
    if (!converged) { warning("optim() failed to converge.  Continuing anyhow.") }

    paramvec(output.config) <- param.optim$par

    # to check things with:
    if (FALSE) {
        true.hts <- hitting.analytic(neighborhoods, G-diag(dG), numcores=numcores )
        range((hts-true.hts)/true.hts,na.rm=TRUE)
        ph <- plot.ht.fn(layer.prefix=file.path(config.dir,config$layer_prefix), nonmissing=nonmissing)
        true.LdL <- params.logistic.setup(paramvec(config),G,update.G,true.hts,zeros,sc.one,layers,transfn,valfn,ndelta,ngamma)
        this.LdL <- params.logistic.setup(paramvec(output.config),G,update.G,hts,zeros,sc.one,layers,transfn,valfn,ndelta,ngamma)
        layout(matrix(1:15,nrow=3))
        plot.nearby( this.LdL$L, paramvec(output.config), fac=.1 )
    }

    if ( converged && all( abs(init.params-paramvec(output.config)) < unlist(config$paramtol) )) { 
        output.config$optim.results <- param.optim
        output.config$optim.nsteps <- this.step
        break
    }
}

if (this.step==config$maxstep) { warning("Maximum number of steps reached.") }

write.json.config( output.config, outfile )
write.full.hts( hts, locs, config$fitted_ht_file )

