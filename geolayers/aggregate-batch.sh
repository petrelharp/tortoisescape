#!/bin/bash
#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:sl230s:ppn=16
#PBS -l walltime=24:00:00
#PBS -l mem=24gb
#PBS -l vmem=120gb
#PBS -l pmem=1500mb

###
# Make the various resolutions we want.

if [ ! -z ${PBS_O_WORKDIR-} ]
then
    source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh
    cd $PBS_O_WORKDIR
else 
    BASEDIR=${1-}
    OUTDIR=${2-}
fi

if [ -z ${BASEDIR-} ] || [ -z ${OUTDIR-} ]
then
    echo "Steps down by factors of 2.  Usage:
    qsub -vBASEDIR=\"expanded/expanded-TIFF\",OUTDIR=\"expanded/\" aggregate-batch.sh
or
    ./aggregate-batch.sh (directory of inputs) (base directory for outputs)
"
    exit
fi

set -eu
set -o pipefail


for RES in 2 4 8 16 32 64 128 256 512
do
    echo "-------------------"
    NEWDIR="${OUTDIR}${RES}x"
    echo "Beginning on ${RES}x and outputting to ${NEWDIR}."
    mkdir -p ${NEWDIR}
    Rscript aggregate-layers.R 2 $NEWDIR $(ls $BASEDIR/*.{grd,tif})
    echo "Done with ${RES}x."
    echo "-------------------"
    BASEDIR=$NEWDIR
done