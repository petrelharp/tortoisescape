#!/usr/bin/Rscript

usage <- "
Create a set of new layers with a different resolution,
  aggregating with means (NAs removed).
Usage:
   Rscript (aggregation factor) (output directory) ( files to aggregate )
"

if (length(commandArgs(TRUE))<3) { stop(usage) }

if (!interactive()) {
    agfact <- as.numeric( commandArgs(TRUE)[1] )
    outdir <- commandArgs(TRUE)[2]
    infiles <- commandArgs(TRUE)[-(1:2)]
} else {
    agfact <- 2.5
    outdir <- "geolayers/TIFF/250x"
    infiles <- list.files("geolayers/TIFF/100x",pattern="*.(tif)|(grd)",full.names=TRUE)
}

require(parallel)
numcores <- detectCores()

require(raster)
rasterOptions(tmpdir=".")

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

mclapply( infiles, function (infile) {
        tryCatch({
                cat("Beginning on ", infile, " .\n")
                raster2crop <- raster(infile)
                outname <- file.path(outdir,gsub("\\.((tif)|(grd))$","",basename(infile)))
                cat("   writing to ", outname, " .\n")
                aggregate(raster2crop,
                            fact=agfact,
                            fun=mean,
                            na.rm=TRUE,
                            filename=outname,
                            overwrite=TRUE)
			 },error=function(e){cat("problem with",infile,"\n")})
        removeTmpFiles()
        cat("Done with ", infile, " .\n")
    }, mc.cores=numcores )

cat("All done.\n")
