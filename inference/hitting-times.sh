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

if [[ -z "${PBS_O_WORKDIR-}" ]]  # not run through pbs
then
    LAYERFILE=${1-}
    PARAMFILE=${2-}
fi

if [[ -z "${LAYERFILE-}" || ! -r "$LAYERFILE" || -z "${PARAMFILE-}" || ! -r "$PARAMFILE" ]]
then
    echo "USAGE:    
        qsub -vLAYERFILE=\"raster-list-file\",PARAMFILE=\"params-file.tsv\" dem-hitting-time.sh 
    or
        ./dem-hitting-time.sh raster-list-file params-file.tsv

    "
    exit 1;
fi

echo "  raster list file:  $LAYERFILE"
echo "  parameter file:  $PARAMFILE"

if [ -e /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh ]
then
    source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh
    cd $PBS_O_WORKDIR
fi

RESLIST="512x 256x 128x 64x 32x 16x" # 8x 4x 2x 1x"

RES="512x"
RESLIST="256x 128x 64x" # 32x 16x 8x 4x" # 2x 1x"

echo "Hitting times at ${RES}."
HTFILE=${RES}/${LAYERFILE}-hitting-times.tsv
# get initial hitting times 
Rscript make-resistance-distances.R ../geolayers/multigrid/${RES}/crm_ ${RES} $LAYERFILE $PARAMFILE analytic ${HTFILE}

for NEXTRES in $RESLIST
do
    echo "-------------------------------------"
    echo "Pushing from ${RES} up to ${NEXTRES}."
    # push these up to the next grid
    Rscript disaggregate-ht.R ../geolayers/multigrid/${RES}/crm_ ../geolayers/multigrid/${NEXTRES}/crm_ ${RES} ${NEXTRES} ${LAYERFILE} ${HTFILE}

    echo "----------------------------"
    echo "Hitting times at ${NEXTRES}."
    HTFILE=${NEXTRES}/${LAYERFILE}-hitting-times.tsv
    Rscript make-resistance-distances.R ../geolayers/multigrid/${NEXTRES}/crm_ ${NEXTRES} $LAYERFILE $PARAMFILE numeric ${HTFILE} ${NEXTRES}/${RES}-${LAYERFILE}-aggregated-hitting-times.tsv 100

    # iterate
    RES=$NEXTRES
done
