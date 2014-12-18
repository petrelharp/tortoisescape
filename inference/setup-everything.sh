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
LAYERDIR="../geolayers/multigrid"

for RES in 512x 256x 128x 64x # 32x 16x 8x 4x 2x 1x
do
    echo $RES
    echo "----------------------------"
    for LAYERS in six-raster-list twelve-raster-list twentyfour-raster-list
    do
        echo $LAYERS
        echo "\n\n----------------------------"
        PREFIX=$(ls ${LAYERDIR}/${RES}/*agp_250.gri | sed -e 's/agp_250.*//')
        echo "   na layer\n\n"
        Rscript make-overlap-na-layer.R ${PREFIX} ${LAYERS}
        echo "\n\n----------------------------"
        echo "   setup G\n\n"
        Rscript setup-real-G.R ${PREFIX} ${RES} ${LAYERS}
        echo \n\n"----------------------------"
        echo "   setup locations\n\n"
        Rscript setup-tort-locs.R ${PREFIX} ${RES} ${LAYERS}
        echo "\n\n----------------------------"
        echo "   save all setup together\n\n"
        Rscript setup-inference.R ${PREFIX} ${RES} ${LAYERS} ${PIFILE}
        echo "\n\n----------------------------"
        echo "Done with $RES : $LAYERS .\n\n"
    done
done
