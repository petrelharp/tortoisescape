# test time complexity
require(Matrix)
source("jacobi.R")
source("resistance-fns.R")

# declare transition/generating matrix
n <- 32
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

######
R <- sparseMatrix(i = c(rightii,upii,leftii,downii),
                j = c(rightjj,upjj,leftjj,downjj),
                x = abs(rnorm(4*n*(n-1),0,1)))
# make the matrix invertible
diags <- rowSums(R) + rexp(nrow(R))/10
D_1 <- sparseMatrix(i = 1:n^2, j = 1:n^2, x = 1/diags)

# pre-calculate stuff  
b <- rnorm(n^2,0,1)
D_R <- (D_1)%*%(R)
D_1b <- (-1)*(D_1)%*%b

# check has the right answer
#  (need to make sure there *is* a right answer first)
solve.x <- solve(R-diag(diags),b)
x0 <- jacobi(D_R,D_1b)
bb <- R %*% x0 - diags*x0
D_1bb <- (-1)*(D_1)%*%bb
x <- jacobi(D_R,D_1bb)
xx <- jacobi(D_R,D_1bb,tol=1e-16)

layout(matrix(1:6,nrow=3))
plot( solve.x, x0 ); abline(0,1)
plot( (R-diag(diags))%*%x0, b ); abline(0,1)
plot( R %*% x - diags*x, bb ); abline(0,1)
plot( bb, ( R %*% x - diags*x - bb ) )
plot( xx, (xx-x) )
plot( bb, ( R %*% xx - diags*xx - bb ) )
# ok, works


###
# crosscheck with hitting.jacobi
n <- 32
loc <- sample.int(n^2,size=1)
b <- rep(-1,n^2)
b[loc] <- 0

R <- sparseMatrix(i = c(rightii,upii,leftii,downii),
                j = c(rightjj,upjj,leftjj,downjj),
                x = abs(rnorm(4*n*(n-1),0,1)))
R <- (R+t(R))/2
diags <- rowSums(R)
D_1 <- sparseMatrix(i = 1:n^2, j = 1:n^2, x = 1/diags)
# set stuff up for hitting calc
hD_1 <- D_1; D_1[,loc] <- 0
hR <- R; hR[,loc] <- 0
hD_R <- hD_1 %*% (hR)
hD_1b <- (-1)*(hD_1)%*%b
hdiags <- diags; hdiags[loc] <- 0

solve.x <- solve( (hR-diag(hdiags))[-loc,-loc], b[-loc] )
solve.x <- c(solve.x[1:(loc-1)],0,solve.x[loc:length(solve.x)])
# jacobi() NOT WORKING; UNCLEAR WHY
jac1.x <- jacobi(hD_R[-loc,-loc],hD_1b[-loc])
jac1.x <- c(jac1.x[1:(loc-1)],0,jac1.x[loc:length(jac1.x)])
jac2.x <- jacobi(hD_R,hD_1b)
hta.x <- hitting.analytic(loc,R-diag(diags))
# hitting.jacobi WORKING; TAKES FOREVER
htj0.x <- hitting.jacobi(loc,R,hts=cbind(solve.x))
htj1.x <- hitting.jacobi(loc,R,hts=cbind(b),kmax=10000)
htj1.x <- hitting.jacobi(loc,R,hts=(htj1.x),kmax=10000)

layout(matrix(1:6,nrow=3))
plot(jac1.x,solve.x); abline(0,1)
plot(jac2.x,solve.x); abline(0,1)
plot(hta.x,solve.x); abline(0,1)
plot(htj0.x,solve.x); abline(0,1)
plot(htj1.x,solve.x); abline(0,1)

layout(matrix(1:4,nrow=2))
plot(((R-diag(diags))%*%solve.x)[-loc])
plot(((R-diag(diags))%*%hta.x)[-loc])
plot(((R-diag(diags))%*%htj0.x)[-loc])
plot(((R-diag(diags))%*%htj1.x)[-loc])

resids <- function (x) { ((hR-diag(hdiags))%*%x - b)[-loc] }
layout(matrix(1:4,nrow=2))
plot(solve.x[-loc],resids(solve.x))
plot(solve.x[-loc],resids(jac1.x))
plot(solve.x[-loc],resids(jac2.x))
plot(solve.x[-loc],resids(hta.x))
