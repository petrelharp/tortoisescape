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

if [[ $NDIRS -lt 1 ]]; then echo "No directories!"; exit; fi

if [[ -e /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh ]]  # run on the cluster
then
    qsub -vBASEDIR=\""$1"\" -t 1-${NDIRS}%32 setup-many-combinations.sh
else
    for SDIR in $(find ${BASEDIR} -mindepth 1 -maxdepth 1 -type d )
    do
        while (( $(jobs 2>&1 | grep -c Running) >= 16 )); do sleep 1; done
        if [ -e ${SDIR}/config.json ]
        then
            Rscript setup-from-json.R ${SDIR}/config.json &
        fi
    done
fi
