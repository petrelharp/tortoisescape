#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:sl230s:ppn=16
#PBS -l walltime=8:00:00
#PBS -l mem=120gb
#PBS -l vmem=120gb
#PBS -l pmem=7500mb

if [[ -z "${PBS_O_WORKDIR-}" ]]  # not run through pbs
then
    SETUP_DIR=${1-}
fi

if [ -e /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh ]
then
    source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh
    cd $PBS_O_WORKDIR
fi

set -eu
set -o pipefail

echo "${SETUP_DIR}"

if [ -z ${SETUP_DIR-} || ! -r ${SETUP_DIR} ]
then
    echo "Usage:  qsub -v'SETUP_DIR=(name of directory)' setup-many-combinations.sh"
fi


for SDIR in $(find ${SETUP_DIR} -type d)
do
    while (( $(jobs 2>&1 | grep -c Running) >= 17 )); do sleep 1; done
    Rscript setup-from-json.R ${SDIR}/config.json ${SDIR}/setup.RData  &
done
