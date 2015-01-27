#!/bin/bash

if [ -e /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh ]
then
    source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh
fi

Rscript -e 'options(width=1000);xx <- t(sapply(commandArgs(TRUE), function (x) { load(x); y <- trust.optim[c("value","iterations","converged")]; y[sapply(y,is.null)]<-NA; unlist(y) } )); xx[order(xx[,1]),]' $(find "$@" -name 'inference-*.RData')
