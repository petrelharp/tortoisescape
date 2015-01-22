##
# Functions to return setup for various optimization problems.
require(trust)

direct.setup <- function (obs.locs, obs.hts, neighborhoods, G, update.G, layers, transfn, valfn, ndelta, ngamma, numcores=getcores()) {
    # here, parameters are:
    #   (T, beta, gamma, delta)
    # where T is an overall shift
    return( function (params) {
            T.shift <- params[1]
            hs <- hitting.sensitivity(params[-1], neighborhoods, G, update.G, layers, transfn, valfn, ndelta, ngamma, do.hessian=TRUE, numcores=numcores )
            hs.grad <- hs$gradient[obs.locs,,]
            new.hts <- (-1)*hs$gradient[obs.locs,,1]  # this works because beta enters as exponential
            resids <- as.vector( new.hts - obs.hts + T.shift )
            dim(hs.grad) <- c( prod(dim(hs.grad)[1:2]), dim(hs.grad)[3] )
            hs.hess <- hs$hessian[obs.locs,,,]
            dim(hs.hess) <- c( prod(dim(hs.hess)[1:2]), prod(dim(hs.hess)[3:4]) )
            ht.hessian <- crossprod( hs.hess, resids )
            dim(ht.hessian) <- c( length(params)-1, length(params)-1 )
            gradient <- 2 * c( sum(resids), crossprod( hs.grad, resids ) )
            hessian <- matrix( 0, nrow=length(params), ncol=length(params) )
            hessian[1,1] <- 2 * length(resids) 
            hessian[-1,1] <- hessian[1,-1] <- 2 * colSums(hs.grad) 
            hessian[-1,-1] <- 2 * ( ht.hessian + crossprod(hs.grad) )
            return( list( value=sum(resids^2), gradient=gradient, hessian=hessian ) )
        } )
}

logistic.trust.setup <- function (G,update.G,hts,zeros,sc.one,layers,transfn,valfn,ndelta,ngamma) {
    # set up for 'trust' region optimization
    # that uses first and second deriv
    #
    # overlaps with params.logistic.setup()
    weightings <- 1/rowMeans(abs(hts),na.rm=TRUE)
    nomitted <- sum( weightings[row(hts)[zeros]] )
    L <- function(params) {
        gamma <- params[1+(1:ngamma)]
        delta <- params[1 + ngamma + (1:ndelta)]
        G@x <- update.G(params)
        dG <- rowSums(G)
        GH <- G %*% hts - dG*hts
        GH[zeros] <- 0
        # function value
        value <- ( sum( weightings*rowSums((GH+sc.one)^2) ) - (nomitted)*sc.one^2 )
        grad <- numeric(length(params))
        hess <- matrix(0,nrow=length(params),ncol=length(params))
        # first, deriv wrt beta:
        grad[1] <- 2 * sum( weightings * rowSums(GH * (GH+sc.one)) )
        hess[1,1] <- grad[1] + 2 * sum( weightings * rowSums(GH^2) )
        # now, gamma:
        ZZ <- (1-transfn(valfn(gamma))) * GH * (
                (1-transfn(valfn(gamma))) * GH +
                (1-2*transfn(valfn(gamma))) * (GH+sc.one) 
            )
        for (kk in 1+(1:ncol(layers))) {
            l.kk <- kk-1
            grad[kk] <- 2 * sum( weightings * rowSums( (layers[,l.kk] * (1-transfn(valfn(gamma))) * GH) * (GH+sc.one)) )
            hess[1,kk] <- hess[kk,1] <- grad[kk] + 2 * sum( weightings * rowSums( (layers[,l.kk] * (1-transfn(valfn(gamma))) * GH^2) ) )
            for (jj in seq(2,length.out=kk-1)) {
                l.jj <- jj-1
                hess[kk,jj] <- hess[jj,kk] <- 2 * ( sum( weightings * rowSums( layers[,l.kk] * layers[,l.jj] * ZZ ) ) )
            }
        }
        # now, delta:
        GLH <- vector(mode="list",length=ncol(layers))  # may be necessary to re-compute each time rather than to store this
        for (kk in 1+ncol(layers)+(1:ncol(layers))) {
            l.kk <- kk-(1+ncol(layers))
            GL <- G
            GL@x <- G@x * ( layers[Gjj,l.kk] + layers[G@i+1L,l.kk] ) * (1-transfn(valfn(delta)[G@i+1L]+valfn(delta)[Gjj]))
            dGL <- rowSums(GL)
            GLH[[kk]] <- GL %*% hts - dGL*hts
            GLH[[kk]][zeros] <- 0
            grad[kk] <- ( 2 * sum( weightings * rowSums( GLH[[kk]] * (GH+sc.one) )  ) )
            hess[1,kk] <- hess[kk,1] <- grad[kk] + 2 * sum( weightings * rowSums( GLH[[kk]] * GH ) )
            for (jj in 1+(1:ncol(layers))) {
                l.jj <- jj-1
                # mixed deriv wrt gamma[jj], delta[kk]
                hess[jj,kk] <- hess[kk,jj] <- 2 * sum( weightings * ( rowSums( layers[,l.jj] * (1-transfn(valfn(gamma))) * GLH[[kk]] * (2*GH+sc.one) ) ) )
            }
            for (ll in seq(2+ncol(layers),length.out=kk-ncol(layers)-1)) {
                # wrt delta[ll] and delta[kk]: product of first partials
                hess[ll,kk] <- hess[kk,ll] <- 2 * sum( weightings * ( rowSums( GLH[[kk]] * GLH[[ll]] ) ) )
            }
            GL@x <- GL@x * 
            for (ll in seq(2+ncol(layers),length.out=kk-ncol(layers)-1)) {
                l.ll <- ll-(1+ncol(layers))
                # wrt delta[ll] and delta[kk]: mixed second deriv
                GL@x <- GL@x * ( layers[Gjj,l.ll] + layers[G@i+1L,l.ll] ) * (1-2*transfn(valfn(delta)[G@i+1L]+valfn(delta)[Gjj]))
                dGL <- rowSums(GL)
                ZZ <- GL %*% hts - dGL * hts
                ZZ[zeros] <- 0
                hess[ll,kk] <- hess[kk,ll] <- hess[ll,kk] + 2 * sum( weightings * rowSums( ZZ * (GH+sc.one) ) )
            }
        }
        return( list( value=value, gradient=grad, hessian=hess ) )
    }
}

update.logistic.trust.setup <- function (L,hts,update.weightings=FALSE) {
    eL <- environment(L)
    assign("hts",hts,eL)
    if (update.weightings) {
        assign("weightings", 1/rowMeans(abs(hts)), eL)
        assign("nomitted", sum( get("weightings",eL)[row(hts)[get("zeros",eL)]] ), eL)
    }
}

params.logistic.setup <- function (init.params,G,update.G,hts,zeros,sc.one,layers,transfn,valfn,ndelta,ngamma) {
    # Given hitting times (hts), return objective function and gradient function of parameters for | G %*% hts + 1 |^2 .
    #
    # setup: evaluating L and dL will change variables they share in a common scope
    # but NOTE: this is actually unnecessary
    #  since L and dL will share the environment created when this function is called.
    # example:
    #  f <- function (x) { e <- environment(); list(u=function (y) { x+y }, v=function(z) {assign("x",z,e)} ) }
    #  gh <- f(3)
    #  gh$u(2)
    #  [1] 5
    #  gh$v(2)
    #  gh$u(2)
    #  [1] 4
    L.env <- new.env()
    assign( "G", G, L.env )
    assign( "update.G", update.G, L.env )
    assign( "hts", hts, L.env )
    assign( "zeros", zeros, L.env )
    assign( "sc.one", sc.one, L.env )
    assign( "transfn", transfn, L.env )
    assign( "valfn", valfn, L.env )
    assign( "ndelta", ndelta, L.env )
    assign( "ngamma", ngamma, L.env )
    assign( "weightings",  1/rowMeans(abs(hts),na.rm=TRUE), L.env )
    assign( "nomitted",  sum( get("weightings",L.env)[row(hts)[zeros]] ), L.env )
    assign("update.aux", function (params,check=TRUE) {
            if ( (!check) || any(params != get("params", L.env ) ) ) { 
                assign("params", params,  L.env )
                evalq( G@x <- update.G(params),  L.env )
                evalq( dG <- rowSums(G),  L.env )
                evalq( GH <- G %*% hts - dG*hts, L.env )
                evalq( GH[zeros] <- 0, L.env )
            }
        }, L.env )

    # initialize
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
        bgrad <- 2 * sum( weightings * rowSums(GH * (GH+sc.one)) )  # wrt beta
        ggrads <- sapply( 1:ncol(layers), function (kk) {  # wrt gamma
                2 * sum( weightings * rowSums( (layers[,kk] * (1-transfn(valfn(gamma))) * GH) * (GH+sc.one)) )
            } )
        dgrads <- sapply( 1:ncol(layers), function (kk) {  # wrt delta
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

###

params.integral.setup <- function () {
    # setup for the integral equation (which does not depend on the functional form, as it does not compute the derivative)
    #   needs to be defined with the things available in the evnironment setup by e.g. params.logistic.setup()
    IL <- function (params) {
        # integral.hts is the mean hitting time of each neighborhood to its boundary,
        #  plus the mean hitting time to each neighborhood (including itself)
        update.aux(params,parent.env(environment()))
        hitting.probs <- get.hitting.probs( G, dG, neighborhoods[nonoverlapping], boundaries[nonoverlapping], numcores=numcores )
        hitting.times <- get.hitting.times( G, dG, neighborhoods[nonoverlapping], boundaries[nonoverlapping], numcores=numcores )
        integral.hts <- do.call( rbind, mclapply( seq_along(neighborhoods[nonoverlapping]), function (k) {
                ihs <- hitting.times[[k]] + hitting.probs[[k]] %*% hts[boundaries[[nonoverlapping[k]]],nonoverlapping] 
                ihs[,k] <- 0
                return(ihs)
            }, mc.cores=numcores ) )
        ans <- sum( ( hts[unlist(neighborhoods[nonoverlapping]),nonoverlapping] / integral.hts - 1 )^2, na.rm=TRUE )
        # ans <- sum( ( hts[unlist(neighborhoods[nonoverlapping]),nonoverlapping] - integral.hts )^2, na.rm=TRUE )
        if (!is.finite(ans)) { browser() }
        return(ans)
    }
    return(IL)
}

###

interp.ht.setup <- function () {
    # For interpolating hitting times to the full grid.
    #
    # Note works on ONE COLUMN at a TIME  (below is ht[locs] not ht[locs,])
    H <- function (ht,locs,zeros,alpha,obs.ht) {
        # ( (G-D)ht + 1 )^T S ( (G-D)ht + 1) + alpha * ( ht[locs] - obs.ht )^T ( ht[locs] - obs.ht )
        # where S = I except S[zeros,zeros]=0
        ht[zeros] <- 0
        z <- G%*%ht - dG*ht + 1
        z[zeros] <- 0
        return( ( sum( z^2 ) )/length(z)  + alpha * sum( ( ht[locs] - obs.ht )^2 ) ) 
    }
    dH <- function (ht,locs,zeros,alpha,obs.ht) {
        # 2 (G-D)^T S ( (G-D) ht + 1 ) + 2 * alpha  * ht[locs]^T ( ht[locs] - obs.ht )[ inserted into big matrix ]
        ht[zeros] <- 0
        z <- G%*%ht - dG*ht + 1
        z[zeros] <- 0
        z <- (crossprod(G,z) - dG*z) / length(z)
        z[zeros] <- 0
        z[locs] <- z[locs] + alpha * crossprod( ht[locs], ht[locs] - obs.ht ) 
        return( 2 * as.vector(z) )
    }
    return( list( H=H, dH=dH ) )
}


d.hts.setup <- function (G,neighborhoods,locs,hts,params,obs.hts) {
    # setup for inferring parameters
    # based instead on minimizing | hts[obs,] - obs.hts |^2 
    R.env <- new.env()
    assign("G", G, R.env)
    assign("neighborhoods", neighborhoods, R.env)
    assign("locs", locs, R.env)
    assign("hts", hts, R.env)
    assign("params", params, R.env)
    assign("obs.hts", obs.hts, R.env)
    assign("update.aux", function (new.params) {
            if (any(new.params!=params)) {
                evalq(G@x <- update.G(params), R.env)
                evalq(hts <- hitting.analytic(neighborhoods,G), R.env)
                evalq(params <- new.params, R.env)
            }
        }, R.env )
    d.hts <- function (new.params) {
        update.aux(new.params)
        return( mean( ( hts[locs,] - obs.hts )^2 ) )
    }
    gr.d.hts <- function (new.params) {
        update.aux(new.params)
        bgrad <- XXX
        return(mean(ans))
    }
}

jacobi.interp.setup <- function ( obs.hts, obs.locs, zeros, alpha=0.2 ) {
    return( function (G,hts,niter) {
        inv.dG <- 1/rowSums(G)
        R <- G
        R@x <- G@x*inv.dG[1L+G@i]
        for (j in 1:niter) {
            hts <- 0.5 * ( hts + inv.dG + R%*%hts )
            hts[obs.locs,] <- (1-alpha) * hts[obs.locs,] + alpha * obs.hts
            hts[zeros] <- 0
        }
        return(hts)
    } )
}

###
# general purpose

gcheck <- function (f,df,params,eps=1e-8) {
    # check gradient: want f1-f0 == df0 if eps is small enough that df0 == df1
    gcheck.fn <- function (k) {
        dirn <- ifelse(seq_along(params)==k,1,0)
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
        results[k,] <- gcheck.fn(k)
        print( results[k,] )
        cat("\n")
    }
    return(results)
}

hcheck <- function (f,params,eps=1e-8,...) {
    # check gradient and hessian: 
    # want f(p+dp)-f(p) = dp * f'(p) [ + (1/2) dp * f''(p) * dp ]
    # and f'(p+dp)-f'(p) = dp * f''(p)
    hcheck.fn <- function (k) {
        dirn <- ifelse(seq_along(params)==k,1,0)
        dp <- eps*dirn
        f0 <- f(params,...)
        f1 <- f(params+dp,...)
        v0 <- f0$value
        v1 <- f1$value
        g0 <- f0$gradient
        g1 <- f1$gradient
        h0 <- f0$hessian
        h1 <- f1$hessian
        c( f0=v0, 
                diff=v1-v0, 
                df0=sum(dp*g0), 
                df1=sum(dp*g1),
                df2.0=sum(dp*g0)+sum(dp*h0%*%dp)/2, 
                df2.1=sum(dp*g1)+sum(dp*h1%*%dp)/2, 
                grad=g0,
                dgrad=g1-g0, 
                ddf0=h0%*%dp, 
                ddf1=h1%*%dp
            )
    }
    results <- sapply( seq_along(params), function (k) {
            # cat("Checking parameter ", k, " : ", names(params)[k], " .\n")
            res <- hcheck.fn(k)
            # print( res )
            # cat("\n")
            return(res)
        } )
    return(results)
}

find.parscale <- function (f, params, parscale, fnscale, eps=0.01, step=10, maxit=3) {
    # guess at an appropriate set of scalings for params,
    # i.e. scalings such that changing each of params by eps*parscale changes f by approximately eps*fnscale
    #  and goes down in one direction
    # or as close as is possible up to changing by 'step'
    fval <- f(params)
    if (missing(fnscale)) { fnscale <- abs(fval) }
    for (k in seq_along(params)) {
        cat(k, " : " )
        for (niter in 1:maxit) {
            dp <- eps*ifelse(seq_along(params)==k,parscale,0)
            fvals <- c( f(params+dp), f(params-dp) )
            log.df <- log10( min( abs( (fvals-fval)/fnscale ), na.rm=TRUE ) ) - log10(eps) 
            # if log.df[k] > 1, we want to make parscale[k] smaller, and vice-versa
            if ( niter==1 ) { 
                dscale <- if ( all( (fvals-fval)/fnscale > 0 ) || log.df > 1 ) { 1/step } else { step }
            }
            if ( ( ( log.df < 0 ) || ( log.df > 1 ) ) && ( log(dscale)*log.df < 0 ) ) { 
                cat(parscale[k], " : ")
                parscale[k] <- parscale[k] * dscale 
            } else {
                cat(parscale[k], " : ", log.df, "\n")
                break
            }
        }
    }
    return(parscale)
}

check.parscale <- function (f, params, parscale, fnscale, eps=0.01) {
    # check if parscale is appropriate,
    # i.e. if near to params, changing each coordinate by an additive factor of parscale changes f by a similar amount
    fval <- f(params)
    if (missing(fnscale)) { fnscale <- abs(fval) }
    results <- sapply( seq_along(params), function (k) {
            dp <- eps*ifelse(seq_along(params)==k,parscale,0)
            fvals <- c( f(params+dp), f(params-dp) )
            log.df <- log10( min( abs( fvals-fval ), na.rm=TRUE ) ) - log10(eps) - log10(fnscale)
        } )
    names(results) <- names(params)
    magnitudes <- floor( log10( results - min(results) + 1 )  ) 
    if ( diff(range(magnitudes)) > 3 ) { warning("Discrepancy in parameter scaling.\n\n") }
    if ( min(magnitudes) < (-1) ) { warning("Some parameter scalings may be too small.\n\n") }
    cat("\nWant these to be all the same order of magnitude, around 1.0:\n\n")
    return(results)
}

compute.plot.nearby <- function (f,params,fac,npoints=20,do.params=seq_along(params)) {
    if (length(fac)==1) { fac <- rep(fac,length(params)) }
    c( 
        lapply( do.params, function (k) {
            parvals <- seq( params[k]-fac[k], params[k]+fac[k], length.out=npoints )
            parmat <- matrix( rep(params,each=npoints), nrow=npoints )
            colnames(parmat) <- names(params)
            parmat[,k] <- parvals
            fvals <- apply( parmat, 1, f )
            return( data.frame(parvals=parvals, fvals=fvals ) )
        } ),
        list(baseval=f(params),params=params,do.params=do.params)
    )
}

plot.nearby <- function (f,params,fac,npoints=20,do.params=seq_along(params),
    computed=compute.plot.nearby(f,params,fac,npoints,do.params), 
    grads=NULL, vlines=NULL, ...) {
    # Make marginal plots of the function f nearby to fac by an additive factor 'fac'
    # makes length(do.params) plots.
    baseval <- computed$baseval
    params <- computed$params
    do.params <- computed$do.params
    for (k in seq_along(do.params)) {
        with( computed[[k]], {
            yrange <- range(fvals,baseval)
            plot( parvals, fvals, ylim=yrange, main=names(params)[do.params[k]], xlim=range(c(parvals,vlines[k])), ... )
            abline(v=c(params[do.params[k]],vlines[k]),col=c("black","red"),lty=c(1,2))
            abline(h=baseval)
            if (length(grads)>0) {
                d.eps <- diff(range(parvals,vlines[k]))/10
                dirn <- (-1)*sign(grads[k])
                arrows( x0=params[do.params[k]], x1=params[do.params[k]]+dirn*d.eps, y0=baseval, y1=baseval+dirn*d.eps*grads[k], col='green' )
            }
        } )
    }
    return(invisible(computed))
}

compute.slice.nearby <- function (f,params,fac,npoints=7,pdir1,pdir2) {
    # take a slice through parameter space in the directions given by pdir1 and pdir2
    # and evaluate f there
    if (length(pdir1)==1) { pdir1 <- ifelse( seq_along(params)==pdir1, 1, 0 ) }
    if (length(pdir2)==1) { pdir2 <- ifelse( seq_along(params)==pdir2, 1, 0 ) }
    res <- matrix(NA,nrow=npoints,ncol=npoints)
    eps <- seq(-fac,fac,length.out=npoints)
    for (j in 1:npoints) { 
        for (k in 1:npoints) {
            res[j,k] <- f( params + eps[j]*pdir1 + eps[k]*pdir2 )
        } 
    }
    return( list( eps=eps, res=res, pdir1=pdir1, params=params ) )
}

plot.slice.nearby <- function (f,params,fac,npoints=7,pdir1,pdir2,col=diverge_hcl(64),computed=compute.slice.nearby(f,params,fac,npoints=7,pdir1,pdir2)) {
    image( computed$eps, computed$eps, computed$res, xlab="pdir1", ylab="pdir2", col=col )
}
