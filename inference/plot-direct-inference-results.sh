#!/bin/bash

if [ $# -lt 1 ]
then
    echo "Usage:  plot-direct-inference-results.sh (name of Rdata file)"
    exit
fi

RMD=$(readlink -f direct-inference-results.Rmd)
DIRNAME=$(dirname $1)
INFILE=$(basename $1)
OUTFILE=$(echo $INFILE|sed -e 's/.[Rr][Dd]ata/.html/')

# note that despite the setwd() below, the .Rmd does *not* behave as if in $DIRNAME.
R -e "setwd(\"${DIRNAME}\");library(knitr);infile=\"$1\";knit2html(\"${RMD}\",output=\"${OUTFILE}\");"

exit
