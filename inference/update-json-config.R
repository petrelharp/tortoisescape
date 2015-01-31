#!/usr/bin/Rscript

usage <- "For each subdirectory in the directory passed in, change the config file so that:
 1. the parameter vector matches the most recent results file
 2. instances of the string 'old_res' are changed to 'new_res' (if present).
Usage:
   Rscript disag-direct-results.R (directory) (old_res) (new_res)
"

if (length(commandArgs(TRUE))<1) { stop(usage) }
basedir <- commandArgs(TRUE)[1]
do.replace <- (length(commandArgs(TRUE))==3)
if (do.replace) {
    old.res <- commandArgs(TRUE)[2]
    new.res <- commandArgs(TRUE)[3]
}

source("input-output-fns.R")
subdirs <- grep( c("figure","cache"), list.dirs(basedir,recursive=FALSE,full.names=TRUE), value=TRUE, inverse=TRUE )
for (subdir in subdirs) {
    config <- read.json.config(file.path(subdir,"config.json"))
    # switch parameters
    result.files <- list.files(subdir,pattern="inference-.*RData",full.names=TRUE)
    last.result <- result.files[ rev(order( file.info(result.files)$mtime )) ][1]
    load(last.result)
    params <- trust.optim$argument
    paramvec(config) <- params
    # gsub resolutions
    if (do.replace) {
        for (k in which(grepl(old.res,config))) {
            config[[k]] <- gsub(old.res,new.res,config[[k]],fixed=TRUE)
        }
    }
    write.json.config(config,file=file.path(subdir,"config.json"))
}
