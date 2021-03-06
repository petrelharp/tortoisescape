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

paramvec <- function (config,vname="params") {
    # access the parameters in a config list as a single vector with nice names
    params <- config[[vname]][c("T","beta","gamma","delta")]
    params$gamma <- c(config[[vname]]$constant_gamma,params$gamma)
    params$delta <- c(config[[vname]]$constant_delta,params$delta)
    names(params$gamma) <- names(params$delta) <- c("constant",config$layer_names)
    return(unlist(params))
}

"paramvec<-" <- function (config,value) {
    # the assignment function corresponding to paramvec()
    value <- as.numeric(value)
    k <- 1
    for (x in c("T","beta","constant_gamma","gamma","constant_delta","delta")) {
        nparams <- length(config$params[[x]]) 
        config$params[[x]][] <- value[ seq(k,length.out=nparams) ]
        k <-  k + nparams
    }
    stopifnot( k == (length(value)+1) )
    return( config )
}

plot.model <- function(params,layer.names,layers,G,update.G,ph,...) {
    # plot the stationary distribution and jump rates for a model
    gamma <- params[2+(1:length(layer.names))]
    stationary.base <- rowSums( layers * gamma[col(layers)] )
    stationary.dist <- ( 1 + exp( -stationary.base ) )  # recall stationary distribution is ** 1/ ** rho( )
    ph( stationary.dist, main="stationary distribution", do.lims=FALSE )
    delta <- params[2+length(layer.names)+(1:length(layer.names))]
    jump.base <- rowSums( layers * delta[col(layers)] )
    G@x <- update.G(params[-1])
    ph( rowSums(G), main="total jump rate", do.lims=FALSE, ... )
}


####
# read/write hitting times etc in standardized way

df.to.sym <- function (x,inds) {
    # takes a three-column data frame x, where the first two have labels in 'inds',
    # and converts it to a matrix with rows and columns in the order of 'inds'
    # and values given by the third column of x,
    # assumed to be symmetric
    ind1 <- match(x[,1],inds)
    ind2 <- match(x[,2],inds)
    mat <- matrix(NA,nrow=length(inds),ncol=length(inds))
    rownames(mat) <- colnames(mat) <- inds
    mat[ cbind(ind1,ind2) ] <- x[,3]
    mat[is.na(mat)] <- t(mat)[is.na(mat)]
    return(mat)
}

sym.to.df <- function ( mat ) {
    # inverse to df.to.sym
    data.frame(
            etort1=rownames(mat)[row(mat)][upper.tri(mat,diag=TRUE)],
            etort2=colnames(mat)[col(mat)][upper.tri(mat,diag=TRUE)],
            dist=mat[upper.tri(mat,diag=TRUE)]
        )
}

read.pairwise.hts <- function ( file, inds=1:n, n=length(inds), upper=TRUE, diag=TRUE ) {
    # read in unstructured hitting times
    # as e.g. output by pwp scripts
    vals <- scan(file) # has UPPER with diagonal
    pimat <- numeric(n^2)
    dim(pimat) <- c(n,n)
    dimnames(pimat) <- list( inds, inds )
    which.tri <- if (upper) { upper.tri } else { lower.tri }
    other.tri <- if (upper) { lower.tri } else { upper.tri }
    pimat[which.tri(pimat,diag=diag)] <- vals
    pimat[other.tri(pimat,diag=FALSE)] <- t(pimat)[other.tri(pimat,diag=FALSE)]
    return(pimat)
}

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

plot.ht.fn <- function (layer.prefix,nonmissing,layer.name="dem_30",
        layer=raster(paste(layer.prefix,layer.name,sep='')),homedir="..",
        sample.loc.file=file.path(homedir,"tort_272_info/geog_coords.RData"),
        default.par.args=list(mar=c(5,4,4,7)+.1),
        default.zlim.fac=1.2) {
    # use this to make a quick plotting function
    require(raster)
    require(rgdal)
    values(layer)[-nonmissing] <- NA # NOTE '-' NOT '!'
    sample.loc.obj <- load(sample.loc.file)
    assign("sample.locs",get(sample.loc.obj))
    sample.locs <- spTransform( sample.locs, CRSobj=CRS(proj4string(layer)))
    orig.locs <- cellFromXY( layer, sample.locs )
    locs <- match(orig.locs,nonmissing)
    ph <- function (x,...,do.lims=TRUE,par.args=default.par.args,zlim.fac=default.zlim.fac) { 
        if (do.lims) {  # restrict to the range observed in observed locations
            lims <- range(x[locs],na.rm=TRUE)
            lims <- mean(lims)+zlim.fac*(lims-mean(lims))
            x <- pmin(lims[2],pmax(lims[1],x))
        }
        values(layer)[nonmissing] <- x
        opar <- par(par.args)  # plotting layers messes up margins
        plot(layer,...)
        points(sample.locs,pch=20,cex=.25)
        par(opar)
    }
    environment(ph) <- new.env()
    assign("sample.locs",sample.locs,environment(ph))
    assign("locs",locs,environment(ph))
    assign("default.par.args",default.par.args,environment(ph))
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

###
# other useful


weighted.median <- function (x,w) { 
    usethese <- w!=0 & !is.na(x) & !is.na(w)
    x <- x[ usethese ]
    w <- w[ usethese ]
    if (length(x)==0) { return( NA ) }
    xord <- order(x)
    x <- x[xord]
    w <- w[xord]
    cw <- cumsum(w)
    k <- min(which(cw>(0.5*cw[length(w)])))
    a <- (0.5*cw[length(w)]-cw[k])/(cw[k+1]-cw[k])
    return( (1-a)*x[k] + a*x[k+1] )
}
