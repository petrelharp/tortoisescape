
setwd("/home/rcf-40/pralph/panfs/tortoises/tortoisescape")

source("resistance-fns.R")
source("profiling-fns.R")
require(raster)

layer.prefix <- c("geolayers/TIFF/masked/")
layer.names <- list.files( layer.prefix, "*.grd" )

system.time( {
    rast <- raster(paste(layer.prefix,layer.names[1],sep=''))
    layer <- values(rast)
    n <- dim(rast)[2]; m <- dim(rast)[1]
} )
#    user  system elapsed 
#  32.531   5.329  37.841 

nvals <- sum(!is.na(layer))

# make smaller
n <- m <- 17500
nvals <- n*m

lsos()
# layer             numeric 4240621208 530077646      NA

rm(layer)
gc()

nlayers <- 4
layers <- vector( mode="list", length=nlayers )
for (k in 1:nlayers) {
    print( system.time( vals <- values( raster(paste(layer.prefix,layer.names[k],sep='')) ) ))  #  32.416   5.275  37.671 
    cat("layer ", k, " nonmissing: ", sum(!is.na(vals)), "\n\n")
    print(system.time( layers[[k]] <- vals[!is.na(vals)][1:(n*m)] ))  #  16.401   2.195  18.597 
}

rm(vals)
gc()

lsos()

# layers  size
#   2   4912554952
#   4   9800000232
#   8   18471436064

print( object.size(layers[[1]]), units="Mb" )
# 2340.5 Mb

system.time( G <- grid.adjacency(n,m,diag=FALSE,symmetric=TRUE) )

print(object.size(G),units="Mb")

#   n     m     user     system   elapsed  object.size
# 13250 13250  804.755   137.067  941.380   4687.7 Mb
# 17500 17500  1410.928  271.167 1681.338   8177.4 Mb

lsos()

#   n     m      G
# 13250 13250 4915433632

layer <- layers[[1]]
system.time( x <- G %*% layer )

#   n     m       user  system elapsed 
# 13250 13250    4.383   0.896   5.276 
# 17500 17500    7.593   1.590   9.183 

system.time( G@x <- G@x * layer[1L+G@i] )
#   n     m         user  system elapsed 
# 17500 17500      7.886   2.178  10.062 

system.time( jj <- rep( seq.int(length(G@p)-1), diff(G@p) ) )
#   n     m         user  system elapsed 
# 17500 17500      18.192   2.463  20.644 

system.time( G@x <- G@x * layer[jj] )
#   n     m         user  system elapsed 
# 17500 17500      6.420   1.585   8.003 

# how many layers can we fit in memory?
128000/2336.5
# 54.78279

# or, how many G matrices?
128000/8177.4
# 15.6529

