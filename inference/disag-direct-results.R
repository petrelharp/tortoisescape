subdirs <- setdiff( list.dirs(".",recursive=FALSE,full.names=FALSE), c("figure","cache") )
for (subdir in subdirs) {
    config <- read.json.config(file.path(subdir,"config.json"))
    result.files <- list.files(subdir,pattern="inference-.*RData",full.names=TRUE)
    last.result <- result.files[ rev(order( file.info(result.files)$mtime )) ][1]
    load(last.result)
    params <- trust.optim$argument
    paramvec(config) <- params
    write.json.config(config,file=file.path(subdir,"config.json"))
}
