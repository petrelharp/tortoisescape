#!/usr/bin/Rscript

# width of square grid
n <- 100
nmap <- matrix(1:n^2,nrow=n)

# random landscape
A <- grid.generator(n)
B <- grid.generator(n)
G <- A+B

# sampling locations
nsamps <- 20
locs.ij <- data.frame( i=sample.int(n,nsamps)-1L, j=sample.int(n,nsamps)-1L )
locs <- ij.to.k(locs.ij,n)

# analytical mean hitting times
true.hts <- hitting.analytic(locs,G)

kk <- sample.int(nsamps,1)
plot( row(nmap), col(nmap), pch=20, col=colorize(true.hts[,kk]) )
points( locs.ij[kk,1], locs.ij[kk,2], pch="*", cex=2 )
points( locs.ij[,1], locs.ij[,2], pch="+", cex=1 )


# Observed hitting times
pairwise.hts <- true.hts[locs,]

# interpolate these
interp.loess <- lapply( seq_along(locs), function (kk) {
            z <- pairwise.hts[,kk]
            loess( z ~ i * j, data=locs.ij )
        } )
all.locs <- data.frame( i=as.vector(row(nmap)), j=as.vector(col(nmap)) )
interp.hts <- sapply( interp.loess, function (interp) {
            predict( interp, newdata=all.locs, span=0.5 )
        } )

# how's this do?
layout(matrix(1:40,nrow=8))
par(mar=c(0,0,0,0))
colorbreaks <- colorize( c(unlist(true.hts),unlist(interp.hts)), return.breaks=TRUE )
invisible( lapply(1:20, function (kk) { 
                    plot( j~i, col=colorize(true.hts[,kk],breaks=colorbreaks), data=all.locs, xlab='', ylab='', xaxt='n', yaxt='n' )
                    plot( j~i, col=colorize(interp.hts[,kk],breaks=colorbreaks), data=all.locs, xlab='', ylab='', xaxt='n', yaxt='n' )
                } ) )


# now try to infer the coefficients:
#  want a,b to minimize 
#     | (a*A+b*B) %*% interp.hts + ones |^2

# these ones are the 'self hitting times' that we don't care about
zeros <- locs + (0:(length(locs)-1))*nrow(interp.hts)

# check this is right
ones <- (-1) * G %*% true.hts
range(ones[-zeros]) # yup
ones[zeros]

# function
f <- function (ab) {
    resids <- ( (ab[1]*A + ab[2]*B) %*% interp.hts + 1 )
    sum( resids^2 - sum(resids[zeros]^2))
}

# gradient
Ath <- crossprod(A,interp.hts)
Athsq <- sum(Ath^2)
Bth <- crossprod(B,interp.hts)
Bthsq <- sum(Bth^2)
f <- function (ab) {
    2 * c( ab[1] * Ath
            ab[1] * crossprod
    c( sum( A  %*% interp.hts , B  %*% interp.hts , 
    resids <- ( (ab[1]*A + ab[2]*B) %*% interp.hts + 1 )
    sum( resids[-zeros]^2 )
}

