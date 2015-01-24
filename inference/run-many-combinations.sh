#!/bin/bash
#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:sl230s:ppn=16
#PBS -l walltime=8:00:00
#PBS -l mem=120gb
#PBS -l vmem=120gb
#PBS -l pmem=7500mb
#PBS -t 1:210%20

if [[ -e /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh ]]
then
    source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh
    cd $PBS_O_WORKDIR
fi

DIRNAME=$(find many-combinations -maxdepth 1 -type 'd' | tail -n +${PBS_ARRAY_INDEX} | head -n 1)

Rscript direct-inference.R ${DIRNAME}/config.json ${DIRNAME}/inference-${PBS_JOBID}_${PBS_ARRAYID}.RData 200

exit
