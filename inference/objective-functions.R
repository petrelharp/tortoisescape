##
# Functions to return setup for various optimization problems.

params.logistic.setup <- function () {
    # This works if the following are defined globally:
    #  G
    #  update.G
    #  hts
    #  weightings
    #  zeros
    #  sc.one
    #  layers
    #  transfn
    #  valfn
    #  ndelta
    #  ngamma

    # setup: evaluating L and dL will change variables they share in a common scope
    L.env <- new.env()
    assign("update.aux", function (params,check=TRUE) {
            if ( (!check) || any(params != get("params", L.env ) ) ) { 
                assign("params", params,  L.env )
                evalq( G@x <- update.G(params),  L.env )
                evalq( dG <- rowSums(G),  L.env )
                evalq( GH <- G %*% hts - dG*hts, L.env )
                evalq( GH[zeros] <- 0, L.env )
            }
        }, L.env )
    # weightings <- ifelse( rowMeans(hts) < quantile(hts,.5), dG, 0 )  # indexes locations; note may overlap with zeros
    # weightings <- ifelse( 1:nrow(hts) %in% locs, 1, 0 )
    assign( "weightings",  1/rowMeans(hts,na.rm=TRUE), L.env )
    assign( "nomitted",  sum( weightings[row(hts)[zeros]] ), L.env )

    evalq( update.aux(init.params,check=FALSE), L.env )

    L <- function (params) {
        update.aux(params)
        ans <- ( sum( weightings*rowSums((GH+sc.one)^2) ) - (nomitted)*sc.one^2 )
        if (!is.finite(ans)) { browser() }
        return(ans)
    }
    dL <- function (params) {
        update.aux(params)
        gamma <- params[1+(1:ngamma)]
        delta <- params[1 + ngamma + (1:ndelta)]
        bgrad <- ( 2 / params[1] )* sum( weightings * rowSums(GH * (GH+sc.one)) )
        ggrads <- sapply( 1:ncol(layers), function (kk) {
                2 * sum( weightings * rowSums( (layers[,kk] * (1-transfn(valfn(gamma))) * GH) * (GH+sc.one)) )
            } )
        dgrads <- sapply( 1:ncol(layers), function (kk) {
                GL <- G
                GL@x <- G@x * ( layers[Gjj,kk] + layers[G@i+1L,kk] ) * (1-transfn(valfn(delta)[G@i+1L]+valfn(delta)[Gjj]))
                dGL <- rowSums(GL)
                GLH <- GL %*% hts - dGL*hts
                GLH[zeros] <- 0
                return( 2 * sum( weightings * rowSums( GLH * (GH+sc.one) )  ) )
            } )
        ans <- ( c(bgrad, ggrads, dgrads) )
        if (any(!is.finite(ans))) { browser() }
        return(ans)
    }
    environment(L) <- environment(dL) <- L.env
    return( list( L=L, dL=dL ) )
}


gcheck <- function (f,df,params,eps=1e-8) {
    # check gradient: want f1-f0 == df0 if eps is small enough that df0 == df1
    gcheck.fn <- function (k) {
        dirn <- ifelse(seq_along(init.params)==k,1,0)
        dp <- eps*dirn
        f0 <- f(params)
        df0 <- df(params)
        f1 <- f(params+dp)
        df1 <- df(params+dp)
        c( f0=f0, difff=f1-f0, df0=sum(dp*df0), df1=sum(dp*df1) )
    }
    results <- matrix(NA, ncol=4, nrow=length(params) )
    for (k in seq_along(params)) {
        cat("Checking parameter ", k, " : ", names(params)[k], " .\n")
        results[k,] <- gcheck(k)
        print( results[k,] )
    }
    return(results)
}

plot.nearby <- function (f,params,fac,npoints=20,...) {
    if (length(fac)==1) { fac <- rep(fac,length(params)) }
    # look at the surface described by f nearby to params
    lapply( seq_along(params), function (k) {
            parvals <- seq( results$par[k]-fac[k], results$par[k]+fac[k], length.out=npoints )
            parmat <- matrix( rep(params,each=npoints), nrow=npoints )
            colnames(parmat) <- names(params)
            parmat[,k] <- parvals
            fvals <- apply( parmat, 1, f )
            yrange <- range(fvals,f(results$par))
            plot( parvals, fvals, ylim=yrange, main=names(params)[k],... )
            abline(v=results$par[k])
            abline(h=f(results$par))
            return( cbind(parvals, fvals=fvals ) )
        } )
}
