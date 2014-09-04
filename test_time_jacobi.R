# test time complexity
require(Matrix)
source("jacobi.R")

nn <- ceiling(exp(log(10)/4 * 1:5))
nreps <- floor(max(nn)/nn)
nn.times <- rep(1,length(nn))

for(i in seq_along(nn)){

  # declare transition/generating matrix
  n <- nn[i]
  rightii <- 1:(n*(n-1))
  rightjj <- rightii + n
  
  upii <- 1:(n^2)
  upii <- upii[-n*(1:n)]
  upjj <- upii + 1
  
  leftii <- (n+1):(n^2)
  leftjj <- leftii - n
  
  downii <- 1:(n^2)
  downii <- downii[-(1+n*(0:(n-1)))]
  downjj <- downii - 1
  
  R <- sparseMatrix(i = c(rightii,upii,leftii,downii),
                    j = c(rightjj,upjj,leftjj,downjj),
                    x = rnorm(4*n*(n-1),0,1))
  D_1 <- sparseMatrix(i = 1:n^2, j = 1:n^2, x = 1/rowSums(abs(R)))

  # free memory
  rm("rightii","rightjj","upii","upjj",
      "leftii","leftjj","downii","downjj")
  gc()
  
  # pre-calculate stuff  
  b = rnorm(n^2,0,1)
  D_R <- (D_1)%*%(-R)
  D_1b <- (D_1)%*%b
  
  # free memory again
  #rm("G_d","D")
  #gc()
  
  t <- system.time(for(k in 1:nreps[i]) jacobi(D_R,D_1b))

  nn.times[i] <- t["user.self"]/nreps[i]
  
}

pdf(file="jacobi-timing.pdf")
plot(nn,nn.times,log='xy',xlab="side length of graph",ylab="seconds to jacobi()")
dev.off()
