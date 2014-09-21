jacobi <- function(D_R,D_1b){

  kmax <- 10000 #maximum iterations
  tol <- 1e-6
  err <- 1 
  x = rnorm(length(b),0,1) # starting vector
  k <- 1

  while(k < kmax && err > tol && err < 2^16){
    x_new <- D_R%*%x + D_1b
    err <- sum((x_new-x)^2)/length(x)
    x <- x_new
    k <- k + 1
  }
  return(x)
}
