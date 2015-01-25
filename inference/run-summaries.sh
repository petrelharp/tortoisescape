#!/bin/bash

if [ -e /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh ]
then
    source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh
fi

Rscript -e 'xx <- t(sapply(commandArgs(TRUE), function (x) { load(x); unlist(trust.optim[c("value","iterations","converged")]) } )); xx[order(xx[,1]),]' $(find "$@" -name 'inference-*.RData')
