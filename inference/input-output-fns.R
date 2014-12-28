##
# configuration 

require(jsonlite)

read.config <- function (...,usage=NULL) {
    # Read in and process the argments passed, and assign them in the parent environment
    # either from the command line or from scan(), depending if interactive or not.
    args <- list(...)
    argvec <- if (interactive()) { cat("Enter command-line parameters ( ", paste(names(args), collapse=", ")," ):\n"); scan(what="character") } else { commandArgs(TRUE) }
    if ( (!is.null(usage)) && ( length(args) != length(argvec) ) ) { stop(usage) }
    for (k in seq_along(args)) { cat( names(args)[k], " : ", argvec[k], "\n" ) }
    argfns <- lapply( args, switch,
                json=read.json.config,
                numeric=as.numeric,
                identity
            )
    argvals <- lapply( seq_along(argfns), function (k) { argfns[[k]](argvec[k]) } )
    names(argvec) <- names(argvals) <- names(args)
    for (k in seq_along(argvals)) {
        assign( names(argvals)[k], argvals[[k]], parent.env(environment()) )
    }
    return( invisible( argvec ) )
}

read.json.config <- function (file) {
    cat("Reading ", file, " .\n")
    con <- openread(file)
    json <- paste(readLines(con, warn = FALSE), collapse = "\n")
    close(con)
    return( fromJSON(json,simplifyMatrix=FALSE) )
}

write.json.config <- function (config,file) {
    json <- toJSON(config,pretty=TRUE)
    cat("Writing to ", file, " .\n")
    con <- openwrite(file)
    writeLines(json,con=con)
    close(con)
}

paramvec <- function (config) {
    # access the parameters in a config list as a single vector with nice names
    params <- config$params[c("beta","gamma","delta")]
    names(params$gamma) <- names(params$delta) <- config$layer_names
    return(unlist(params))
}

"paramvec<-" <- function (config,value) {
    # the assignment function corresponding to paramvec()
    value <- as.numeric(value)
    k <- 1
    for (x in c("beta","gamma","delta")) {
        nparams <- length(config$params[[x]]) 
        config$params[[x]][] <- value[ seq(k,length.out=nparams) ]
        k <-  k + nparams
    }
    stopifnot( k == (length(value)+1) )
    return( config )
}


####
# read/write hitting times etc in standardized way

write.full.hts <- function ( hts, locs, file ) {
    # write out hitting times as a matrix
    colnames( hts ) <- locs
    cat("Writing output to ", file, " .\n")
    write.table( hts, file=file, row.names=FALSE )
}

read.full.hts <- function (file,locs) {
    # file should be in the following format:
    #   "144" "23" "201"
    #    23.0 45.1 120.2
    #    33.3 31.3 110.8
    #    42.0 85.8 124.2
    # and will check the first row corresponds to 'locs'
    hts.locs <- as.numeric( scan(file,nlines=1,what="char") )
    stopifnot( all( hts.locs == locs ) )
    matrix( scan( file, skip=1 ), ncol=length(hts.locs), byrow=TRUE )
}

write.sub.hts <- function ( hts, file ) {
    hts.df <- data.frame(
            row=as.vector(row(hts)), # [upper.tri(hts,diag=TRUE)],
            col=as.vector(col(hts)), # [upper.tri(hts,diag=TRUE)],
            DISTANCE=as.vector(hts)  # [upper.tri(hts,diag=TRUE)]
        )
    cat("Writing out to ", file, " .\n")
    write.table( noisy.df, file=file, row.names=FALSE )
}

read.sub.hts <- function ( file, locs ) {
    # Read in the three-column format of hitting times with row and column labels
    #  TO-DO: should be using tortoise IDs, not row and column numbers like currently.
    ht.df <- read.table( file, header=TRUE, stringsAsFactors=FALSE)
    ht <- matrix( NA, nrow=length(locs), ncol=length(locs) )
    ht[ cbind( ht.df$row, ht.df$col ) ] <- ht.df$DISTANCE
    return( ht )
}

# miscellaneous utility functions

selfname <- function (x) { names(x) <- x; return(x) }

# load a file, but not into the global environment, rather, into a list.
load.to.list <- function (file) { e <- environment(); n <- load(file,envir=e); names(n) <- n; lapply( n, get, envir=e ) }

# use to open stdin/stdout or process substitution things correctly
#   from  http://stackoverflow.com/questions/15784373/process-substitution
openread <- function(arg) {
    if (arg %in% c("-", "/dev/stdin","stdin")) {
       stdin()
    } else if (grepl("^/dev/fd/", arg)) {
       fifo(arg, open = "r")
    } else {
       file(arg, open = "r")
    }
}
openwrite <- function(arg) {
    if (arg %in% c("-", "/dev/stdout","stdout")) {
       stdout()
    } else if (grepl("^/dev/fd/", arg)) {
       fifo(arg, open = "w")
    } else {
       file(arg, open = "w")
    }
}

###
# use for debugging noninteractive stuff
# options(error=print.and.dump)

print.and.dump <- function () {
 cat(paste("Error in \"", paste(commandArgs(),collapse=' '), "\": dumping frames.\n")); dump.frames(to.file = TRUE); q(status=1)
} 

###
# random seed

getset.seed <- function () {
    new.seed <- as.integer(runif(1)*2e9)
    set.seed(new.seed)
    cat("Setting seed to ", new.seed, "\n")
    return(new.seed)
}

##
# plotting whatnot

plot.ht.fn <- function (layer.prefix,layer.name="dem_30",nonmissing,layer=raster(paste(layer.prefix,layer.name,sep='')),homedir="..",par.args=list(mar=c(5,4,4,7)+.1)) {
    # use this to make a quick plotting function
    require(raster)
    values(layer)[-nonmissing] <- NA # NOTE '-' NOT '!'
    load(paste(homedir,"tort.coords.rasterGCS.Robj",sep='/'))
    ph <- function (x,...) { 
        values(layer)[nonmissing] <- x
        opar <- par(par.args)  # plotting layers messes up margins
        plot(layer,...)
        points(tort.coords.rasterGCS,pch=20,cex=.25)
        par(opar)
    }
    environment(ph) <- new.env()
    assign("tort.coords.rasterGCS",tort.coords.rasterGCS,environment(ph))
    return(ph)
}

colorize <- function (x, nc=32, colfn=function (n) rainbow_hcl(n,c=100,l=50), zero=FALSE, trim=0, breaks, return.breaks=FALSE) {
    if (is.numeric(x) & trim>0) {
        x[ x<quantile(x,trim,na.rm=TRUE) ] <- quantile(x,trim,na.rm=TRUE)
        x[ x>quantile(x,1-trim,na.rm=TRUE) ] <- quantile(x,1-trim,na.rm=TRUE)
    }
    if (missing(breaks) & is.numeric(x)) {
        if (zero) {
            breaks <- seq( (-1)*max(abs(x),na.rm=TRUE), max(abs(x),na.rm=TRUE), length.out=nc )
        } else {
            breaks <- seq( min(x,na.rm=TRUE), max(x,na.rm=TRUE), length.out=nc )
        }
        x <- cut(x,breaks=breaks,include.lowest=TRUE)
    } else {
        x <- factor(x)
    }
    if (return.breaks) {
        return(breaks)
    } else {
        return( colfn(nlevels(x))[as.numeric(x)] )
    }
}

