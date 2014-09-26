#!/bin/bash

# do all the setup

for RES in 500x 100x 10x masked
do
    for LAYERS in six-raster-list twelve-raster-list twentyfour-raster-list
    do
        PREFIX=$(ls ../geolayers/TIFF/${RES}/*agp_250.gri | sed -e 's/agp_250.*//')
        Rscript make-overlap-na-layer.R ${PREFIX} ${LAYERS}
        Rscript setup-real-G.R ${PREFIX} ${LAYERS}
        Rscript setup-tort-locs.R ${PREFIX}
    done
done

for RES in 500x 100x
do
    for LAYERS in six-raster-list twelve-raster-list twentyfour-raster-list
    do
        Rscript setup-inference.R ${PREFIX} ${RES} ${LAYERS}
    done
done
