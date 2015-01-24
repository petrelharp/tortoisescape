#!/bin/bash
#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:sl230s:ppn=16
#PBS -l walltime=48:00:00
#PBS -t 1-210%20
#PBS -l mem=24gb
#PBS -l vmem=24gb
#PBS -l pmem=1500mb
## max on sl230s:
# #PBS -l mem=120gb
# #PBS -l vmem=120gb
# #PBS -l pmem=7500mb
# # -t denotes the job array limits

###
# notes on memory usage:
#   256x needs < 6GB memory

if [[ -e /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh ]]
then
    source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh
    echo "Changing to ${PBS_O_WORKDIR}"
    cd $PBS_O_WORKDIR
fi

pwd

JOBID=$(echo $PBS_JOBID | sed -e 's/[^0-9].*//')
DIRNAME=$(find many-combinations -mindepth 1 -maxdepth 1 -type 'd' | tail -n +${PBS_ARRAYID} | head -n 1)  
CONFIG="${DIRNAME}/config.json"
OUTPUT="${DIRNAME}/inference-${JOBID}_${PBS_ARRAYID}.RData"
echo "JOBID: ${JOBID}"
echo "PBS_ARRAYID: ${PBS_ARRAYID}"
echo "DIRNAME: ${DIRNAME}"
echo "CONFIG: ${CONFIG}"
echo "OUTPUT: ${OUTPUT}"

Rscript direct-inference.R ${CONFIG} ${OUTPUT} 200

exit
