#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:dl165:ppn=24
#PBS -l walltime=200:00:00
#PBS -l mem=48000mb
#PBS -l vmem=48000mb
#PBS -l pmem=2000mb

# First, do
#  qsub -vARGS="dem-layer-list" setup-multigrid.sh

set -eu
set -o pipefail

if [ -e /home/rcf-40/pralph/cmb/lib/openblas-usc-R-setup.sh ]
then
    source /home/rcf-40/pralph/cmb/lib/openblas-usc-R-setup.sh
    cd $PBS_O_WORKDIR
fi

LAYERLIST="dem-layer-list"
PARAMFILE="params-dem-layer-list.tsv"

RES="512x"
RESLIST="256x 128x" # 64x 32x 16x 8x 4x" # 2x 1x"

echo "Hitting times at ${RES}."
# get initial hitting times 
Rscript make-resistance-distances.R ../geolayers/multigrid/${RES}/crm_ ${RES} $LAYERLIST $PARAMFILE analytic

for NEXTRES in $RESLIST
do
    echo "-------------------------------------"
    echo "Pushing from ${RES} up to ${NEXTRES}."
    # push these up to the next grid
    Rscript disaggregate-ht.R ../geolayers/multigrid/${RES}/crm_ ../geolayers/multigrid/${NEXTRES}/crm_ ${RES} ${NEXTRES} ${LAYERLIST} ${RES}/${LAYERLIST}-hitting-times.tsv 2

    echo "----------------------------"
    echo "Hitting times at ${NEXTRES}."
    Rscript make-resistance-distances.R ../geolayers/multigrid/${NEXTRES}/crm_ ${NEXTRES} $LAYERLIST $PARAMFILE numeric ${NEXTRES}/${RES}-${LAYERLIST}-aggregated-hitting-times.tsv 100

    # iterate
    RES=$NEXTRES
done
