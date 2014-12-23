#!/bin/bash

# Make some random parameter values,
# compute hitting times,
# then fit logistic model from them.

BASEPARAMS="_test_six_layers/six-params.tsv"

# Create random parameter values
R --vanilla --slave << EOF
nsims <- 2
sim.dirs <- file.path("_test_six_layers",paste("test",formatC(1e6*runif(nsims),width=6,format="d",flag="0"),sep="_"))
param.df <- read.table("$BASEPARAMS",header=TRUE)
for (k in 1:nsims) {
    dir.create(sim.dirs[k],showWarnings=FALSE,recursive=TRUE)
    new.params <- 2*rbeta( ncol(param.df)-2, 40, 40 ) - 1
    param.df[,-(1:2)] <- new.params[col(param.df[,-(1:2)])]
    outfile <- file.path(sim.dirs[k],"six-params.tsv")
    cat("Writing to ", outfile, " .\n")
    write.table( param.df, file=outfile, row.names=FALSE, quote=FALSE )
}
EOF

for SIMDIR in "_test_six_layers/test_*"
do
    mkdir -p ${SIMDIR}/256x
    HTFILE=${SIMDIR}/256x/six-raster-list-hitting-times-full.tsv
    Rscript make-resistance-distances.R ../geolayers/multigrid/256x/crm_ 256x six-raster-list ${SIMDIR}/six-params.tsv analytic ${HTFILE}
    TRUTH_START=_test_six_layers/inferred-six-params-true-start.tsv
    Rscript fit-logistic-model.R ../geolayers/multigrid/256x/crm_ 256x six-raster-list ${HTFILE} ${SIMDIR}/six-params.tsv ${TRUTH_START} &
    RANDOM_START=_test_six_layers/inferred-six-params-random-start.tsv
    Rscript fit-logistic-model.R ../geolayers/multigrid/256x/crm_ 256x six-raster-list ${HTFILE} ${BASEPARAMS} ${RANDOM_START} &
    wait
done
