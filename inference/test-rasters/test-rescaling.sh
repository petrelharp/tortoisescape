#!/usr/bin/Rscript

set -eu
set -o pipefail

for RES in 500x 400x 300x 200x 100x
do
    PREFIX=${RES}/${RES}_
    Rscript ../make-overlap-na-layer.R ${PREFIX} test-layers
    Rscript ../setup-real-G.R ${PREFIX} test-layers ${RES}
    (cd ..; Rscript setup-tort-locs.R test-rasters/${PREFIX} test-rasters/${RES})
done

for RES in 500x 400x 300x 200x 100x
do
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
layout(t(1:3))
for (k in 1:ncol(new.hts.500x)) {
    par(mar=c(5,4,4,8)+.1)
    ph(new.hts.500x[,k])
    ph(new.hts.400x[,k])
    ph(new.hts.400x[,k]-new.hts.500x[,k])
    if (is.null(locator(1))) { break }
}
EOF

echo $RCODE 
