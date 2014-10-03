#!/usr/bin/Rscript

set -eu
set -o pipefail

for x in grd gri
do
    cp ../../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_dem_30.$x  100x/100x_dem_30.$x
done

Rscript -e "require(raster); dem <- raster('100x/100x_dem_30'); for (fact in 2:5) { aggregate(dem,fact=fact,fun=mean,na.rm=TRUE,file=paste(fact*100,'x/',fact*100,'x_dem_30',sep=''),overwrite=TRUE) }"

for RES in 500x 400x 300x 200x 100x
do
    PREFIX=${RES}/${RES}_
    Rscript ../make-overlap-na-layer.R ${PREFIX} test-layers
    Rscript ../setup-real-G.R ${PREFIX} test-layers ${RES}
    (cd ..; Rscript setup-tort-locs.R test-rasters/${PREFIX} test-rasters/${RES})
done

for RES in 500x 400x 300x
do
    PREFIX=${RES}/${RES}_
    Rscript ../make-resistance-distances.R ${PREFIX} ${RES} test-layers test-params.tsv analytic
done

Rscript ../disaggregate-ht.R 500x/500x_ 100x/100x_ 500x 100x 500x/test-layers-hitting-times.tsv 5
Rscript ../disaggregate-ht.R 400x/400x_ 100x/100x_ 400x 100x 400x/test-layers-hitting-times.tsv 4
Rscript ../disaggregate-ht.R 300x/300x_ 100x/100x_ 300x 100x 300x/test-layers-hitting-times.tsv 3

read -r -d '' RCODE << EOF
source("resistance-fns.R")
new.hts.500x <- read.table("100x/500x-aggregated-hitting-times.tsv",header=TRUE)
new.hts.400x <- read.table("100x/400x-aggregated-hitting-times.tsv",header=TRUE)
new.hts.300x <- read.table("100x/300x-aggregated-hitting-times.tsv",header=TRUE)
load("100x/100x_nonmissing.RData")
ph <- plot.ht.fn("100x/100x_","dem_30",nonmissing,homedir='../../') 
pdf(file="comparisons.pdf",width=10,height=8,pointsize=10)
layout(matrix(1:6,nrow=2,byrow=TRUE))
for (k in 1:ncol(new.hts.500x)) {
    par(mar=c(5,4,4,8)+.1)
    ph(new.hts.500x[,k],main="500x")
    ph(new.hts.400x[,k],main="400x")
    ph(new.hts.300x[,k],main="300x")
    ph(new.hts.400x[,k]-new.hts.500x[,k],main="400x-500x")
    ph(new.hts.300x[,k]-new.hts.500x[,k],main="300x-500x")
    ph(new.hts.300x[,k]-new.hts.400x[,k],main="300x-400x")
    if (interactive()) if (is.null(locator(1))) { break }
}
dev.off()
EOF

Rscript <(echo $RCODE)
