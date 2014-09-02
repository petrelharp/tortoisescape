require(gsl)
# see http://www.bibsonomy.org/bibtex/2e3ede022c042c956308ca4055b7af1ef/peter.ralph
#  


# angular basis:
# angular(x,k)[z] =
#  if k odd  :  cos( arg(x) * (k+1)/2 )
#  if k even :  sin( arg(x) * k/2 )
# with origin at x0
angular <- function (x,mm,x0=0) {
    if (is.null(dim(x))) { dim(x) <- c(1,length(x)) }
    if (any(x0!=0)) { x <- sweep(x,2,x0,"-") }
    arg <- atan2( x[,1], x[,2] )
    if (mm==0) {
        rep.int(1.0,NROW(x))
    } else if ((mm %% 2) == 1) { # odd
        cos( arg * (mm+1)/2 )
    } else {
        sin( arg * mm/2 )
    }
}

# radial part:
#   Bessel functions
#  see http://en.wikipedia.org/wiki/Fourier%E2%80%93Bessel_series
#   and http://en.wikipedia.org/wiki/Hankel_transform
# use nu=1

radial <- function (x,k,m,x0=c(0,0),a=1) {
    if (is.null(dim(x))) { dim(x) <- c(1,length(x)) }
    if (any(x0!=0)) { x <- sweep(x,2,x0,"-") }
    rad <- rowSums( x^2 )
    bessel_Jnu(nu=m,x=rad * bessel_zero_Jnu(nu=m,k)/a )
}

# "fb" for Fourier-Bessel
fbfun <- function (x,k,mm,m=ceiling(mm/2),x0=c(0,0),a=1) {
    if (is.null(dim(x))) { dim(x) <- c(1,length(x)) }
    if (any(x0!=0)) { x <- sweep(x,2,x0,"-") }
    angular(x,mm=mm) * radial(x,k=k,m=m,a=a)
}

fbbasis <- function (x,mmax,x0=c(0,0),a=1) {
    # inefficent
    if (is.null(dim(x))) { dim(x) <- c(1,length(x)) }
    if (any(x0!=0)) { x <- sweep(x,2,x0,"-") }
    fbb <- do.call( cbind, lapply( 0:(2*mmax), function (mm) {
            do.call( cbind, lapply( 1:ceiling(mmax/2), function (k) {
                        fbfun( x, k=k, mm=mm, a=a )
                    } ) )
            } ) )
    sweep( fbb, 2, sqrt(colSums(fbb^2)), "/" )
}
