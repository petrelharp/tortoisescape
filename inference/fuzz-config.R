#!/usr/bin/Rscript

usage <- "Fuzz the parameters in a config file, creating another one.
Usage:
    Rscript fuzz-config.R (input config) (fuzz factor) (output config)

"

source("input-output-fns.R")

read.config( config="json", fuzz.fac="numeric", outfile="character", usage=usage )

params <- paramvec(config)
cat("Old parameters:\n    ")
cat(paste(names(params),params,sep=" : ", collapse=", "), "\n")
params <- params * exp( rnorm(length(params)) * fuzz.fac )
cat("New parameters:\n    ")
cat(paste(names(params),params,sep=" : ", collapse=", "), "\n")
paramvec(config) <- params
write.json.config(config, outfile)
