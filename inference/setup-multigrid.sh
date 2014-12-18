#!/bin/bash
#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:sl230s:ppn=16
#PBS -l walltime=200:00:00
#PBS -l mem=120gb
#PBS -l vmem=120gb
#PBS -l pmem=7500mb

if [[ -z "${PBS_O_WORKDIR-}" ]]  # not run through pbs
then
    LAYERFILE=$1
fi

if [[ -z "${LAYERFILE-}" || ! -r "$LAYERFILE" ]]
then
    echo "USAGE:    \
        qsub -vLAYERFILE=\"raster-list-file\" setup-multigrid.sh
    or\
        ./setup-multigrid.sh raster-list-file
    "
    exit 1;
fi

echo "raster list file:  $LAYERFILE"

if [ -e /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh ]
then
    source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh
    cd $PBS_O_WORKDIR
fi

set -eu
set -o pipefail

RESLIST="512x 256x 128x 64x 32x 16x 8x 4x" # 2x 1x"
PIFILE="../pairwisePi/alleleCounts_1millionloci.pwp"

for RES in $RESLIST
do
    echo "$RES,  $LAYERFILE"
    mkdir -p $RES
    echo "----------------------------"
    PREFIX=$(ls ../geolayers/multigrid/${RES}/*agp_250.gri | sed -e 's/agp_250.*//')
    echo "   na layer"
    echo "----------------------------"
    Rscript make-overlap-na-layer.R ${PREFIX} ${LAYERFILE}
    echo "   setup G"
    echo "----------------------------"
    Rscript setup-real-G.R ${PREFIX} ${RES} ${LAYERFILE}
    echo "   setup locations"
    echo "----------------------------"
    Rscript setup-tort-locs.R ${PREFIX} ${RES} ${LAYERFILE}
done

for RES in $RESLIST
do
    echo "$RES,  $LAYERFILE"
    echo "----------------------------"
    PREFIX=$(ls ../geolayers/multigrid/${RES}/*agp_250.gri | sed -e 's/agp_250.*//')
    Rscript setup-inference.R ${PREFIX} ${RES} ${LAYERFILE} ${PIFILE}
done
