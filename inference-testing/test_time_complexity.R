# test time complexity
require(Matrix)

nn <- ceiling(exp(log(10)/4 * 1:13))
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
  
  G <- sparseMatrix(i = c(rightii,upii,leftii,downii,1:(n^2)),
                    j = c(rightjj,upjj,leftjj,downjj,1:(n^2)),
                    x = rnorm(5*n^2,0,1))
  # free memory
  rm("rightii","rightjj","upii","upjj",
      "leftii","leftjj","downii","downjj")
  gc()
  
  b = rnorm(n^2,0,1)
  
  t <- system.time(for(k in 1:nreps[i]) solve(G,b))

  nn.times[i] <- t["user.self"]/nreps[i]
  
}

pdf(file="sparse-LU-timing.pdf")
plot(nn,nn.times,log='xy',xlab="side length of graph",ylab="seconds to solve()")
dev.off()
