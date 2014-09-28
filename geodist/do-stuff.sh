#!/bin/bash
# setup for doing landscape layer analysis

for x in $(cat ../inference/twentyfour-raster-list); do Rscript get-environmental-distance.R ../geolayers/TIFF/500x/500x_ 500x $x max; done &
for x in $(cat ../inference/twentyfour-raster-list); do Rscript get-environmental-distance.R ../geolayers/TIFF/500x/500x_ 500x $x mean; done &
