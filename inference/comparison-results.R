#!usr/bin/Rscript

usage <- "Use a fitted model to make predictions in a way comparable to others.
Usage:
    Rscript comparison-results.R (name of results .RData file) (name of summary directory) [name of output .RData file] [name of output .png file]
and the outputs default to the input file, with 'inference-' replaced with 'comparison-' and appended with name of summary directory.
"

argvec <- if (interactive()) { scan(what='char') } else { commandArgs(TRUE) }
if (length(argvec)<2) { stop(usage) }

source("resistance-fns.R")
require(parallel)
numcores<-getcores()
require(raster)

infile <- argvec[1]
summary.dir <- argvec[2]
outfile <- if (length(argvec)>2) { argvec[3] } else {
    gsub( "[.]R[dD]ata$", paste("-",gsub("/","_",summary.dir),".RData",sep=''), 
        gsub("inference-", "comparison-", infile) ) }
outpng <- if (length(argvec)>3) { argvec[4] } else { gsub("[.]R[dD]ata$", ".png", outfile) }

if (!file.exists(infile)){ stop(paste(infile," doesn't exist.")) }
load(infile)

# setup for using *everyone*, on habitat only
config.file <- file.path("summaries",summary.dir,"config.json")
if (!file.exists(config.file)){ stop(paste(config.file," doesn't exist.")) }
config <- read.json.config(config.file)
for (x in config$setup_files) { 
    setup.file <- file.path(dirname(config.file),x)
    if (!file.exists(setup.file)){ stop(paste(setup.file," doesn't exist.")) }
    load(setup.file)
}
config$reference_inds <- row.names( sample.locs )

# and for this model
if (!file.exists(trust.optim$config.file)) { stop(paste(trust.optim$config.file)," doesn't exist.") }
local.config <- read.json.config( trust.optim$config.file )
layer.names <- local.config$layer_names
paramvec(local.config) <- trust.optim$argument
# get centering and scaling used for the layers (HACK)
local.env <- new.env()
for (x in local.config$setup_files) { 
    setup.file <- file.path(dirname(trust.optim$config.file),x)
    if (!file.exists(setup.file)){ stop(paste(setup.file," doesn't exist.")) }
    load(setup.file,envir=local.env) 
}
layer.center <- get("layer.center",local.env)
layer.scale <- get("layer.scale",local.env)

# get locally-specified layers in globally-specified way
use.files <- layer.files[match(layer.names,layer.file.names)]
names(use.files) <- layer.names

## find if other layers have additional NA values, and block these later
new.nalayer <- nalayer
for (lf in use.files) {
    other <- raster(lf) 
    new.nalayer <- mask( new.nalayer, other )
}
block.these <- match( which( is.na(values(new.nalayer)) & ! is.na(values(nalayer)) ), nonmissing )

# load up the layers
if (length(layer.names)>0) { 
    layers <- sapply( layer.names, function (ln) {
                values( raster( paste(full.layer.prefix,ln,sep='') ) )[nonmissing]
            } )
} else {
    layers <- matrix(0,nrow=nrow(G),ncol=0)
}
stopifnot(nrow(layers)==nrow(G))
# ADD the constant layer
layers <- cbind( 1, layers )
layer.names <- c( "constant", layer.names )
ndelta <- ngamma <- length(layer.names)
layers <- sweep(sweep(layers,2,layer.center,"-"),2,layer.scale,"/")

G@x <- update.G(paramvec(local.config)[-1])

hts <- hitting.analytic( neighborhoods, G, numcores=numcores, blocked=block.these )
hts <- hts[locs,]
rownames(hts) <- colnames(hts) <- row.names(sample.locs)

trust.optim.results <- trust.optim[ "converged" ]
    
save( infile, config.file, config, local.config, trust.optim.results, hts, pimat, file=outfile )

## and, make a nice plot
pcs <- read.csv(file.path(dirname(config.file),dirname(config$divergence_file),"pcs.csv"),header=TRUE)
pc.cols <- adjustcolor( ifelse( pcs$PC1[match(rownames(pimat),pcs$etort)][row(pimat)] < 0, "purple", "blue" ), 0.5 )
# remove duplicates: these are (etort-156 / etort-296 ) and (etort-143 / etort-297)
dup.inds <- match( c( "etort-296", "etort-297" ), rownames(pimat) )
omit.comparisons <- ( pcs$PC1[match(rownames(pimat),pcs$etort)][row(pimat)] * pcs$PC1[match(colnames(pimat),pcs$etort)][col(pimat)] < 0 )
omit.comparisons <- ( omit.comparisons | (row(pimat) %in% dup.inds) | (col(pimat) %in% dup.inds) )
# and omit self comparisons
omit.comparisons <- ( omit.comparisons | (row(pimat) == col(pimat))  )

fitted.asym <- paramvec(local.config)[1] + hts
fitted <- paramvec(local.config)[1] + (hts+t(hts))/2 
resids <- (fitted - pimat)
resids[omit.comparisons] <- NA
fitted[omit.comparisons] <- NA
fitted.asym[omit.comparisons] <- NA

# weight residuals by 1 / number of other samples within 25km
geodist.tab <- read.csv(file.path(dirname(config.file),dirname(config$divergence_file),"geog_distance.csv"),header=TRUE,stringsAsFactors=FALSE)
geodist.tab <- subset( ( geodist.tab, etort1 %in% rownames(pimat) ) & ( geodist.tab, etort1 %in% colnames(pimat) ) )
gg <- cbind( match(geodist.tab$etort1,rownames(pimat)), match(geodist.tab$etort2,rownames(pimat)) )
geodist <- pimat
geodist[] <- NA
geodist[ gg ] <- geodist.tab[,3]
geodist[is.na(geodist)] <- t(geodist)[is.na(geodist)]
nearby.weights <- 1 / rowSums( geodist$distance < 25e3 )

med.resids <- apply(resids,1,weighted.median,w=nearby.weights)

dem <- raster( file.path( dirname( gsub("/inference/.*","/inference",getwd()) ), "visualization", "dem_30.gri" ) )
dem.locs <- spTransform( sample.locs, CRSobj=CRS(proj4string(dem)) )

png(file=outpng, width=10*288, height=5*288, res=288, pointsize=10)
layout(t(1:2))
plot( fitted, pimat, pch=20, cex=0.5, col=pc.cols, ylim=range(pimat[!omit.comparisons]), 
   xlab="fitted hitting times", ylab="observed divergence" )
abline(0,1)
plot( dem, main="median residuals", legend=FALSE )
contour( dem, add=TRUE )
points( dem.locs, pch=21, cex=sqrt(abs(med.resids)/5e3), col=ifelse(med.resids>0,"red","blue"), bg=adjustcolor(ifelse(med.resids>0,"red","blue"),.5) )
lvals <- (c(2,1,-1,-2)*5e3)
legend("topleft", pch=21, pt.cex=sqrt(abs(lvals/5e3)), col=ifelse((lvals/5e3)>0,"red","blue"), pt.bg=adjustcolor(ifelse((lvals/5e3)>0,"red","blue"),.5),
    legend=sprintf("%0.1e",lvals) )
dev.off()
