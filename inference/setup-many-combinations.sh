#!/bin/bash

if [[ $# -lt 1 ]]
then
    echo "Usage:  \
        ./setup-many-combinations.sh (name of directory)\
"
    exit
fi

BASEDIR=$1
NDIRS=$(find $1 -mindepth 1 -maxdepth 1 -type 'd' | wc -l)
echo "Parsing ${NDIRS} directories in ${BASEDIR}"

if [[ -e /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh ]]  # run on the cluster
then
    qsub -vBASEDIR=\""$1"\" -t 1-${NDIRS}%32 setup-many-combinations.sh
else
    for SDIR in $(find ${BASEDIR} -type d -mindepth 1 -maxdepth 1)
    do
        while (( $(jobs 2>&1 | grep -c Running) >= 16 )); do sleep 1; done
        Rscript setup-from-json.R ${SDIR}/config.json ${SDIR}/setup.RData  &
    done
fi
