#!/bin/bash

if [ $# -lt 1 ]
then
    echo "Usage:  plot-direct-inference-results.sh (name of Rdata file)"
    exit
fi

OUTFILE=$(echo $1|sed -e 's/.[Rr][Dd]ata/.html/')

R -e "library(knitr);infile=\""$1"\";knit2html('direct-inference-results.Rmd',output=\""${OUTFILE}"\");"

exit
