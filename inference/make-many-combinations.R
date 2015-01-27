layer.info <- read.csv("../geolayers/raster-effects-from-lit.csv",header=TRUE,stringsAsFactors=FALSE)
layer.info$twentyfour <- ( layer.info$X24.rasters != "" )
layer.info$twelve <- ( layer.info$X12.rasters != "" )

layers <- subset(layer.info,twelve)

# always use these
always.use <- c( "dem_30", "agp_250" )
# additional layers
n.additional <- 4
# combinations
lcombns <- layers$layer.name[combn(setdiff(1:nrow(layers),match(always.use,layers$layer.name)),n.additional)]
dim(lcombns) <- c( n.additional, length(lcombns)/n.additional )
layer.combns <- cbind( sapply( always.use, rep, ncol(lcombns) ), t(lcombns) )
stopifnot( all( max( apply( layer.combns, 1, table ) ) == 1 ) )

combn.names <- apply( layer.combns, 1, paste, collapse="-" )

json.config.fn <- function (layer.names) {
  return( paste('
{
    "description" : "Configuration for direct inference on 256x grid.",
    "layer_names" : [ "', paste( layer.names, collapse='", "' ),  '" ],
    "layer_prefix" : "../../../geolayers/multigrid/256x/crm_",
    "mask.layer": [ "mask_crew_dem_2K_sea" ],
    "setup_files" : [ "setup.RData" ],
    "sample_locs" : [ "../../../tort_180_info/tort.coords.rasterGCS.Robj" ],
    "divergence_file" : [ "../../../tort_180_info/alleleCounts_1millionloci.pwp.csv" ],
    "params" : {
            "T"              : [ 3.5e5 ],
            "beta"           : [ 2.0 ],
            "constant_gamma" : [ 0.0 ],
            "gamma"          : [ ', paste( c("-1.0", "1.0", rep("0.0",length(layer.names)-2)), collapse=', '), ' ],
            "constant_delta" : [ 0.0 ],
            "delta"          : [ ', paste( c("-1.0", "1.0", rep("0.0",length(layer.names)-2)), collapse=', '), ' ]
        },
    "param_scale" : {
            "T"              : [ 3.5e5 ],
            "beta"           : [ 1.0 ],
            "constant_gamma" : [ 1.0 ],
            "gamma"          : [ ', paste( rep("1.0",length(layer.names)), collapse=", " ), ' ],
            "constant_delta" : [ 1.0 ],
            "delta"          : [ ', paste( rep("1.0",length(layer.names)), collapse=", " ), ' ]
        }
}
', sep='' ) )
}

for (k in seq_along(combn.names)) {
    outdir <- paste("many-combinations/",combn.names[k],sep='')
    dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
    cat( json.config.fn(layer.combns[k,]), file=file.path(outdir,"config.json") )
}

