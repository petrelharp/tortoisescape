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

qsub -vBASEDIR=\""$1"\" -t 1-${NDIRS}%32 setup-many-combinations.sh
