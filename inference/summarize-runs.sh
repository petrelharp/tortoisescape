#!/bin/bash

if [ -e /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh ]
then
    source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh
fi

Rscript <(echo 'options(width=1000);
    xx <- t(sapply(commandArgs(TRUE), function (x) { 
        load(x); 
        y <- trust.optim[c("value","iterations","converged")];
        y[sapply(y,is.null)]<-NA;
        y[[1]] <- sqrt(y[[1]])/(nrow(trust.optim$pimat)*(ncol(trust.optim$pimat)-1))
        unlist(y) } )); 
    xx[order(xx[,1]),]') \
    $(find "$@" -name 'inference-*.RData')
