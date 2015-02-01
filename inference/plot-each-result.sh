#!/bin/bash

if [ $# -lt 1 ]
then
    echo "Usage: ./plot-each-result.sh (directory to look for results in)"
    exit
fi

for x in ${1}/*/inference-*RData; do y=$(echo $x|sed -e 's/RData/html/'); if [ ! -e $y ]; then ./plot-direct-inference-results.sh $x; fi; done
