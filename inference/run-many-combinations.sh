#!/bin/bash

if [[ $# -lt 2 ]] || [[ ! -d $1 ]]
then
    echo "Usage:\
        ./run-many-combinations.sh (name of directory) (number of iterations)\
"
    exit
fi

NDIRS=$(find $1 -mindepth 1 -maxdepth 1 -type 'd' | wc -l)

"qsub -vBASEDIR="$1",MAXIT=$2 -t 1-${NDIRS}%20 run-many-combinations.pbs"
qsub -vBASEDIR="$1",MAXIT=$2 -t 1-${NDIRS}%20 run-many-combinations.pbs
