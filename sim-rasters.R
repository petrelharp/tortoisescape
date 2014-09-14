#!/usr/bin/Rscript

# Simulate up some rasters to test one

require(raster)

dims <- c(1e3,2e3)
nlayers <- 3

xlist <- replicate( nlayers, { x <- raster(nrow=dims[1],ncol=dims[2]); values(x) <- rexp(prod(dims)); x }, simplify=FALSE )

lapply(seq_along(xlist), function (k) writeRaster(xlist[[k]], file=paste("test-raster-",k,sep='') ) )

xbrick <- brick( replicate( nlayers, { x <- raster(nrow=dims[1],ncol=dims[2]); values(x) <- rexp(prod(dims)); x }, simplify=FALSE ) )

writeRaster(xbrick,file="test-raster-brick")
