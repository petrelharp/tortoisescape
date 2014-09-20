wt.cov <- function (x) {
    nr <- (ncol(x)/2)
    vals <- x[,2*(1:nr)-1]
    wts <- x[,2*(1:nr)]
    wt.prods <- matrix(NA,nrow=nr,ncol=nr)
    wt.denoms <- matrix(NA,nrow=nr,ncol=nr)
    wt.sums <- matrix(NA,nrow=nr,ncol=nr)
    for (k in 1:nr) {
        for (j in 1:k) {
            wt.prods[k,j] <- sum( vals[,k]*vals[,j]*sqrt(wts[,k]*wts[,j]), na.rm=TRUE )
            wt.denoms[k,j] <- sum( sqrt(wts[,k]*wts[,j]), na.rm=TRUE )
            wt.sums[k,j] <- sum( vals[,k] * sqrt(wts[,k]*wts[,j]), na.rm=TRUE )
            wt.sums[j,k] <- sum( vals[,j] * sqrt(wts[,k]*wts[,j]), na.rm=TRUE )
        }
    }
    wt.prods/wt.denoms - (wt.sums/wt.denoms)*(t(wt.sums)/wt.denoms)
}

sub.cov <- function (x) {
    nr <- (ncol(x)/2)
    vals <- ( runif(nrow(x)) > x[,2*(1:nr)-1]/x[,2*(1:nr)] )
    wts <- (x[,2*(1:nr)]>0)
    wt.prods <- matrix(NA,nrow=nr,ncol=nr)
    wt.denoms <- matrix(NA,nrow=nr,ncol=nr)
    wt.sums <- matrix(NA,nrow=nr,ncol=nr)
    for (k in 1:nr) {
        for (j in 1:k) {
            wt.prods[k,j] <- sum( vals[,k]*vals[,j]*sqrt(wts[,k]*wts[,j]), na.rm=TRUE )
            wt.denoms[k,j] <- sum( sqrt(wts[,k]*wts[,j]), na.rm=TRUE )
            wt.sums[k,j] <- sum( vals[,k] * sqrt(wts[,k]*wts[,j]), na.rm=TRUE )
            wt.sums[j,k] <- sum( vals[,j] * sqrt(wts[,k]*wts[,j]), na.rm=TRUE )
        }
    }
    wt.prods/wt.denoms - (wt.sums/wt.denoms)*(t(wt.sums)/wt.denoms)
}

x <- read.table("test.data",header=TRUE)
print( wt.cov(x) )

print(sub.cov(x))

read.table(pipe("cat test.data | awk -f subsampled-covariance.awk"),header=FALSE)

y <- read.table("test.binary.data",header=FALSE)
print( wt.cov(y) )

read.table(pipe("cat test.binary.data | awk -f subsampled-covariance.awk"),header=FALSE)

