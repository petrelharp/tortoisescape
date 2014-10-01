#!/bin/bash

# get initial hitting times on 500x grid
#   makes 500x/six-raster-list-hitting-times.tsv
Rscript make-resistance-distances.R ../geolayers/TIFF/500x/500x_ 500x ../inference/six-raster-list simple-init-params-six-raster-list.tsv analytic

# push these up to 100x grid
Rscript disaggregate-ht.R ../geolayers/TIFF/500x/500x_ ../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_ 500x 100x 500x/six-raster-list-hitting-times.tsv

# now use those to find hitting times on 100x grid
Rscript make-resistance-distances.R ../geolayers/TIFF/100x/crop_resampled_masked_aggregated_100x_ 100x ../inference/six-raster-list simple-init-params-six-raster-list.tsv CG 500x/six-raster-list-hitting-times.tsv

