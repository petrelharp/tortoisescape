###
# Modified from pracma to allow sparse matrices


gmres <- function (A, b, x0 = rep(0, length(b)), errtol = 1e-06, kmax = length(b) + 
    1, reorth = 1) 
{
    b <- as.matrix(b)
    n <- length(b)
    if (nrow(A) != n || ncol(A) != n) 
        stop("Matrix 'A' must be square and compatible with 'b'.")
    h <- zeros(kmax)
    v <- zeros(n, kmax)
    c <- zeros(kmax + 1, 1)
    s <- zeros(kmax + 1, 1)
    normF <- function(x) norm(as.matrix(x), type = "F")
    x <- as.matrix(x0)
    if (norm(x, "F") != 0) {
        r <- b - A %*% x
    }
    else {
        r <- b
    }
    rho <- norm(r, "F")
    g <- rho * eye(kmax + 1, 1)
    errtol <- errtol * norm(b, "F")
    error <- c()
    error <- c(error, rho)
    niter <- 0
    if (rho < errtol) 
        return(list(x = x, error = error, niter = niter))
    v[, 1] <- r/rho
    beta <- rho
    k <- 0
    while (rho > errtol && k < kmax) {
        k <- k + 1
        v[, k + 1] <- as.vector( A %*% v[, k] )
        normav <- normF(v[, k + 1])
        for (j in 1:k) {
            h[j, k] <- t(v[, j]) %*% v[, k + 1]
            v[, k + 1] <- v[, k + 1] - h[j, k] * v[, j]
        }
        h[k + 1, k] <- normF(v[, k + 1])
        normav2 <- h[k + 1, k]
        if ((reorth == 1 && normav + 0.001 * normav2 == normav) || 
            reorth == 3) {
            for (j in 1:k) {
                hr <- t(v[, j]) %*% v[, k + 1]
                h[j, k] <- h[j, k] + hr
                v[, k + 1] = v[, k + 1] - hr * v[, j]
            }
            h[k + 1, k] <- normF(v[, k + 1])
        }
        if (h[k + 1, k] != 0) 
            v[, k + 1] <- v[, k + 1]/h[k + 1, k]
        if (k > 1) 
            h[1:k, k] <- .givapp(c[1:(k - 1)], s[1:(k - 1)], h[1:k, k], k - 1)
        nu <- normF(h[k:(k + 1), k])
        if (nu != 0) {
            c[k] <- Conj(h[k, k]/nu)
            s[k] <- -h[k + 1, k]/nu
            h[k, k] <- c[k] * h[k, k] - s[k] * h[k + 1, k]
            h[k + 1, k] <- 0
            g[k:(k + 1)] <- .givapp(c[k], s[k], g[k:(k + 1)], 1)
        }
        rho <- abs(g[k + 1])
        error <- c(error, rho)
    }
    y <- qr.solve(h[1:k, 1:k], g[1:k])
    x <- x0 + v[1:n, 1:k] %*% y
    return(list(x = x, error = error, niter = k))
}

zeros <- function (n, m = n) 
{
    # this should NOT return a Matrix (hella slow!!)
    stopifnot(is.numeric(n), length(n) == 1, is.numeric(m), length(m) == 1)
    n <- floor(n)
    m <- floor(m)
    if (n <= 0 || m <= 0) 
        return(matrix(0, 0, 0))
    else return(matrix(0, n, m))
}

eye <- function (n, m = n) 
{
    stopifnot(is.numeric(n), length(n) == 1, is.numeric(m), length(m) == 1)
    n <- floor(n)
    m <- floor(m)
    nm <- max(n,m)
    if (n <= 0 || m <= 0) 
        return(matrix(NA, 0, 0))
    else return(base::diag(1,n,m))
}

.givapp <- function (c, s, v_in, k) 
{
    v_rot <- v_in
    for (i in 1:k) {
        w1 <- c[i] * v_rot[i] - s[i] * v_rot[i + 1]
        w2 <- s[i] * v_rot[i] + Conj(c[i]) * v_rot[i + 1]
        v_rot[i:(i + 1)] <- c(w1, w2)
    }
    v_rot
}
