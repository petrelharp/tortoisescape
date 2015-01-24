#!/bin/bash

if [[ $# -lt 1 ]] || [[ ! -d $1 ]]
then
    echo "Usage:\
        ./run-many-combinations.sh (name of directory)\
"
    exit
fi

NDIRS=$(find $1 -mindepth 1 -maxdepth 1 -type 'd' | wc -l)

qsub -vBASEDIR="$1" -t 1-${NDIRS}%20 run-many-combinations.pbs
