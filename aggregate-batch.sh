#!/bin/bash
#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:sl230s:ppn=16
#PBS -l walltime=200:00:00
#PBS -l mem=120gb
#PBS -l vmem=120gb
#PBS -l pmem=7500mb

###
# Make the various resolutions we want.

if [ -e /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh ]
then
    source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh
    cd $PBS_O_WORKDIR
fi

set -eu
set -o pipefail

BASEDIR="geolayers/masked"

for RES in 2 4 8 16 32 64 128 256 512
do
    echo "-------------------"
    echo "Beginning on ${RES}x."
    NEWDIR="geolayers/multigrid/${RES}x"
    mkdir -p ${NEWDIR}
    Rscript aggregate-layers.R 2 $NEWDIR $(ls $BASEDIR/*.{grd,tif})
    echo "Done with ${RES}x."
    echo "-------------------"
    BASEDIR=$NEWDIR
done
