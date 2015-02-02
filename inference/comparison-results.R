#!usr/bin/Rscript

usage <- "Use a fitted model to make predictions in a way comparable to others.
Usage:
    Rscript comparison-results.R (name of results .RData file) [name of output .RData file] [name of output .png file]
and the outputs default to the input file, but with -comparison.RData and -comparison.png appended.
"

argvec <- if (interactive()) { scan(what='char') } else { commandArgs(TRUE) }
if (length(argvec)<1) { stop(usage) }

source("resistance-fns.R")
require(parallel)
numcores<-getcores()
require(raster)
require(rgdal)

infile <- argvec[1]
load(infile)
outfile <- if (length(argvec)>1) { argvec[2] } else { gsub("[.]R[dD]ata$", "-comparison.RData", infile) }
outpng <- if (length(argvec)>2) { argvec[3] } else { gsub("[.]R[dD]ata$", ".png", outfile) }

# setup for using *everyone*, on habitat only
config.file <- "summaries/habitat-only/config.json"
config <- read.json.config(config.file)
for (x in config$setup_files) { load(file.path(dirname(config.file),x)) }
config$reference_inds <- row.names( sample.locs )

# and for this model
local.config <- read.json.config( trust.optim$config.file )
layer.names <- local.config$layer_names
paramvec(local.config) <- trust.optim$argument

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
                scale( values( raster( paste(full.layer.prefix,ln,sep='') ) )[nonmissing] )
            } )
} else {
    layers <- matrix(0,nrow=nrow(G),ncol=0)
}
stopifnot(nrow(layers)==nrow(G))
# ADD the constant layer
layers <- cbind( 1, layers )
layer.names <- c( "constant", layer.names )
ndelta <- ngamma <- length(layer.names)

G@x <- update.G(paramvec(local.config)[-1])

hts <- hitting.analytic( neighborhoods, G, numcores=numcores, blocked=block.these )
hts <- hts[locs,]
rownames(hts) <- colnames(hts) <- row.names(sample.locs)

trust.optim.results <- trust.optim[ "converged" ]
    
save( infile, config.file, config, local.config, trust.optim.results, hts, pimat, file=outfile )

## and, make a nice plot
pcs <- read.csv(file.path(dirname(config.file),dirname(config$divergence_file),"pcs.csv"),header=TRUE)
omit.comparisons <- ( pcs$PC1[match(rownames(pimat),pcs$etort)][row(pimat)] * pcs$PC1[match(colnames(pimat),pcs$etort)][col(pimat)] < 0 )
pc.cols <- adjustcolor( ifelse( pcs$PC1[match(rownames(pimat),pcs$etort)][row(pimat)] < 0, "purple", "blue" ), 0.5 )
# remove duplicates: these are (etort-156 / etort-296 ) and (etort-143 / etort-297)
dup.inds <- match( c( "etort-296", "etort-297" ), rownames(pimat) )
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
geodist <- pointDistance(sample.locs,lonlat=FALSE)
nearby.weights <- 1 / rowSums( geodist < 25e3 )

med.resids <- apply(resids,1,weighted.median,w=nearby.weights)

dem <- raster( file.path( dirname( gsub("/inference/.*","/inference",getwd()) ), "visualization", "dem_30.gri" ) )
dem.locs <- spTransform( sample.locs, CRSobj=CRS(proj4string(dem)) )

png(file=outpng, width=4*288, height=4*288, res=288, pointsize=10)
layout(t(1:2))
plot( fitted, pimat, pch=20, cex=0.5, col=pc.cols, ylim=range(pimat[!omit.comparisons]), 
   xlab="fitted hitting times", ylab="observed divergence" )
abline(0,1)
plot( dem, main="median residuals" )
contour( dem, add=TRUE )
points( dem.locs, pch=21, cex=sqrt(abs(med.resids)/5e3), col=ifelse(med.resids>0,"red","blue"), bg=adjustcolor(ifelse(med.resids>0,"red","blue"),.5) )
lvals <- (c(2,1,-1,-2)*5e3)
legend("topleft", pch=21, pt.cex=sqrt(abs(lvals/5e3)), col=ifelse((lvals/5e3)>0,"red","blue"), pt.bg=adjustcolor(ifelse((lvals/5e3)>0,"red","blue"),.5),
    legend=sprintf("%0.1e",lvals) )
dev.off()
