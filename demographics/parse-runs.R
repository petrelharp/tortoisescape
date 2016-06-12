# Get the results from all the iterations
library(jsonlite)

all.iter.dirs <- grep( "run_.*iter_", list.dirs( "." ), value=TRUE )
have.score <- file.exists( file.path(all.iter.dirs,"model.score") )
iter.dirs <- all.iter.dirs[have.score]
all.params.list <- lapply( iter.dirs, function (itd) { fromJSON(file.path(itd,"params.json")) } )
all.params <- do.call( rbind, lapply( all.params.list, function (x) {
                    dim(x$refugia.coords) <- c(1,4)
                    as.data.frame(x)
        } ) )
all.params$score <- sapply( file.path(iter.dirs,"model.score"), scan, quiet=TRUE )
all.params <- all.params[ order(all.params$score,decreasing=FALSE), ]

pairs(all.params)
