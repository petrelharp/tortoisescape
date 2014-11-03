#!/bin/bash
#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:sl230s:ppn=16
#PBS -l walltime=200:00:00
#PBS -l mem=120gb
#PBS -l vmem=120gb
#PBS -l pmem=7500mb

# do all the setup

if [ -e /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh ]
then
    source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh
    cd $PBS_O_WORKDIR
fi

set -eu
set -o pipefail

PIFILE="../pairwisePi/alleleCounts_1millionloci.pwp"

for RES in 500x 100x 10x # masked
do
    echo $RES
    echo "----------------------------"
    for LAYERS in six-raster-list twelve-raster-list twentyfour-raster-list
    do
        echo $LAYERS
        echo "----------------------------"
        PREFIX=$(ls ../geolayers/TIFF/${RES}/*agp_250.gri | sed -e 's/agp_250.*//')
        echo "   na layer"
        echo "----------------------------"
        Rscript make-overlap-na-layer.R ${PREFIX} ${LAYERS}
        echo "   setup G"
        echo "----------------------------"
        Rscript setup-real-G.R ${PREFIX} ${RES} ${LAYERS}
        echo "   setup locations"
        echo "----------------------------"
        Rscript setup-tort-locs.R ${PREFIX} ${RES} ${LAYERS}
    done
done

for RES in 500x 100x
do
    echo $RES
    echo "----------------------------"
    for LAYERS in six-raster-list twelve-raster-list twentyfour-raster-list
    do
        echo $LAYERS
        echo "----------------------------"
        PREFIX=$(ls ../geolayers/TIFF/${RES}/*agp_250.gri | sed -e 's/agp_250.*//')
        Rscript setup-inference.R ${PREFIX} ${RES} ${LAYERS} ${PIFILE}
    done
done
