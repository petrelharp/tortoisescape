#!/bin/bash
#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:sl230s:ppn=16
#PBS -l walltime=200:00:00
#PBS -l mem=120gb
#PBS -l vmem=120gb
#PBS -l pmem=7500mb

# do all the setup

set -eu
set -o pipefail

if [ -e /home/rcf-40/pralph/cmb/lib/openblas-usc-R-setup.sh ]
then
    source /home/rcf-40/pralph/cmb/lib/openblas-usc-R-setup.sh
    cd $PBS_O_WORKDIR
fi

for RES in 500x 100x 10x # masked
do
    for LAYERS in six-raster-list twelve-raster-list twentyfour-raster-list
    do
        PREFIX=$(ls ../geolayers/TIFF/${RES}/*agp_250.gri | sed -e 's/agp_250.*//')
        Rscript make-overlap-na-layer.R ${PREFIX} ${LAYERS}
        Rscript setup-real-G.R ${PREFIX} ${LAYERS} ${RES}
        Rscript setup-tort-locs.R ${PREFIX} ${RES}
    done
done

for RES in 500x 100x
do
    for LAYERS in six-raster-list twelve-raster-list twentyfour-raster-list
    do
        PREFIX=$(ls ../geolayers/TIFF/${RES}/*agp_250.gri | sed -e 's/agp_250.*//')
        Rscript setup-inference.R ${PREFIX} ${RES} ${LAYERS}
    done
done
