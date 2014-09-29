#!/bin/bash


for x in $(cat ../inference/twentyfour-raster-list)
do
    qsub -vARGS="get-environmental-distance.R ../geolayers/TIFF/10x/crop_resampled_masked_aggregated_10x_ 10x $x mean" ../cluster/single-skeleton.pbs
    qsub -vARGS="get-environmental-distance.R ../geolayers/TIFF/10x/crop_resampled_masked_aggregated_10x_ 10x $x max" ../cluster/single-skeleton-dl165.pbs
done

exit

