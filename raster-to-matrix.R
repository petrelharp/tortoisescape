#!/usr/bin/R

args <- commandArgs(TRUE)

if (length(args) < 1 || length(args) %% 2 != 0) {
    cat("Usage:\n Rscript raster-to-matrix.R base-constant [ name-of-raster-file coefficient-for-raster ]*\n")
}

if (interactive()) {
    args <- c("0", "test-raster-1.grd", "1", "test-raster-2.grd", "-1.4", "test-raster-3.grd", "2.3")
}

# transform from values to migration rates
transform <- function (x) { 1/(1+exp((-1)*x)) }

const <- as.numeric(args[1])
args <- args[-1]
dim(args) <- c(2,length(args)/2)
raster.files <- args[1,]
if (any(!file.exists(raster.files))) { stop("Can't read raster file.") }
alpha <- as.numeric(args[2,])
nrasters <- length(raster.files)

# combine layers
rast <- raster(raster.files[1]) * alpha[1] + const
if (nrasters>1) for (k in 2:nrasters) {
    rast <- rast + alpha[k] * raster(raster.files[k])
}
rast <- transform(rast)

# put this into a generator matrix
G <- grid.adjacency(nrow(rast),ncol(rast),diag=FALSE,symmetric=TRUE)

# 


