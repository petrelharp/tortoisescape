#!/bin/bash

if [[ $# -lt 3 ]] || [[ ! -d $1 ]]
then
    echo "Usage:\
        ./run-many-combinations.sh (name of directory) (number of iterations) (number of concurrent jobs)\
"
    exit
fi

NDIRS=$(find $1 -mindepth 1 -maxdepth 1 -type 'd' | wc -l)

echo "qsub -vBASEDIR="$1",MAXIT=$2 -t 1-${NDIRS}%$3 run-many-combinations.pbs"
qsub -vBASEDIR="$1",MAXIT=$2 -t 1-${NDIRS}%$3 run-many-combinations.pbs
