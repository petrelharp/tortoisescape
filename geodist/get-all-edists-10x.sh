#!/bin/bash


for x in $(cat ../inference/twentyfour-raster-list)
do
    qsub -vARGS="get-environmental-distance.R ../geolayers/TIFF/10x/crop_resampled_masked_aggregated_10x_ 10x $x" ../cluster/single-skeleton.pbs
done

exit

