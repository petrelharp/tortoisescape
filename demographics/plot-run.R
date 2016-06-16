# (re)-plot previously-run runs
library(jsonlite)
library(sp)
library(rgeos)
library(ape)
library(raster)
library(msarg)
source("grid-refugia-fns.R")

iter.dirs <- if ( interactive() ) { scan(what='char') } else { commandArgs(TRUE) }

for (outdir in iter.dirs) {

    basedir <- dirname(outdir)
    base.params <- fromJSON( file.path(basedir,"params.json") )
    if (is.null(base.params$hab.fact)) { base.params$hab.fact <- 32 }
    params <- fromJSON( file.path(outdir,"params.json") )
    dist.df <- read.csv( file.path( basedir, "dist-df.csv" ), header=TRUE, stringsAsFactors=FALSE )
    sim.dist <- scan( file.path(outdir,"sim-distances.csv") )
    # mean-square difference:
    model.score <- scan( file.path(outdir,"model.score") )
    tree.output <- trees_from_ms( file.path(outdir,"msoutput.txt") )


    # make the plots
    # pdf(file=file.path(outdir,"trees-and-things.pdf"),width=10,height=5,pointsize=10)
    png(file=file.path(outdir, "trees-and-things-%02d.png"), 
        width=10*144, height=5*144, pointsize=10, res=144)

    plot_everything( base.params, params, dist.df, sim.dist, model.score, tree.output, label=outdir, plot.ntrees=4 )

    dev.off()
}
