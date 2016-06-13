# Get the results from all the iterations
library(jsonlite)

# function to examine output of a single run
get_results <- function (basedir=".") {
    all.iter.dirs <- grep( "run_.*iter_", list.dirs( basedir ), value=TRUE )
    have.score <- file.exists( file.path(all.iter.dirs,"model.score") )
    iter.dirs <- all.iter.dirs[have.score]
    all.params.list <- lapply( iter.dirs, function (itd) { fromJSON(file.path(itd,"params.json")) } )
    all.params <- do.call( rbind, lapply( all.params.list, function (x) {
                        dim(x$refugia.coords) <- c(1,4)
                        if (is.null(x$hab.fact)) { x$hab.fact <- 32 }
                        as.data.frame(x)
            } ) )
    all.params$score <- sapply( file.path(iter.dirs,"model.score"), scan, quiet=TRUE )
    all.params$dir <- iter.dirs
    all.params <- all.params[ order(all.params$score,decreasing=FALSE), ]
    return(all.params)
}

all.params <- get_results()
last.all.params <- if (file.exists("all-results.csv")) { read.csv("all-results.csv", header=TRUE) } else { NULL }
all.params <- merge( all.params, last.all.params, all.x=TRUE, all.y=TRUE )
all.params <- all.params[ order(all.params$score,decreasing=FALSE), ]

write.csv(all.params, file="all-results.csv", row.names=FALSE)

# pairs(all.params)

# look at variation between runs with the same parameter

# ntrees=500 for these -- looks good.
same <- get_results("run_all_same")
sd(same$score) #= 163433
sd(same$score)/mean(same$score)  # = .037
