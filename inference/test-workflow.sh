#!/bin/bash
#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:sl230s:ppn=16
#PBS -l walltime=200:00:00
#PBS -l mem=120gb
#PBS -l vmem=120gb
#PBS -l pmem=7500mb

if [[ ! -z "${PBS_O_WORKDIR-}" ]]  # run through pbs
then
    source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh
    cd $PBS_O_WORKDIR
fi

# Make some random parameter values,
# compute hitting times,
# then fit logistic model from them.

### Prototype:
if [[ '' ]]
then
    Rscript make-resistance-distances.R ../geolayers/multigrid/256x/crm_ 256x six-raster-list test_six_layers/six-params.tsv analytic test_six_layers/256x/six-raster-list-hitting-times-full.tsv
    Rscript fuzz-hitting-times.R ../geolayers/multigrid/256x/crm_ 256x  six-raster-list test_six_layers/256x/six-raster-list-hitting-times-full.tsv 0.00 test_six_layers/256x/six-raster-list-sim-0_00-hts.tsv
    Rscript fuzz-hitting-times.R ../geolayers/multigrid/256x/crm_ 256x  six-raster-list test_six_layers/256x/six-raster-list-hitting-times-full.tsv 0.005 test_six_layers/256x/six-raster-list-sim-0_005-hts.tsv
    Rscript fuzz-hitting-times.R ../geolayers/multigrid/256x/crm_ 256x  six-raster-list test_six_layers/256x/six-raster-list-hitting-times-full.tsv 0.01 test_six_layers/256x/six-raster-list-sim-0_01-hts.tsv
    Rscript fuzz-hitting-times.R ../geolayers/multigrid/256x/crm_ 256x  six-raster-list test_six_layers/256x/six-raster-list-hitting-times-full.tsv 0.05 test_six_layers/256x/six-raster-list-sim-0_05-hts.tsv
fi

BASEDIR="test_six_layers"
BASEPARAMS="${BASEDIR}/six-params.tsv"
NSIMS=2

# Create random parameter values
R --vanilla --slave << EOF
nsims <- ${NSIMS}
sim.dirs <- file.path("${BASEDIR}",paste("test",formatC(1e6*runif(nsims),width=6,format="d",flag="0"),sep="_"))
param.df <- read.table("${BASEPARAMS}",header=TRUE)
for (k in 1:nsims) {
    dir.create(sim.dirs[k],showWarnings=FALSE,recursive=TRUE)
    new.params <- 2*rbeta( ncol(param.df)-2, 40, 40 ) - 1
    param.df[,-(1:2)] <- new.params[col(param.df[,-(1:2)])]
    outfile <- file.path(sim.dirs[k],"six-params.tsv")
    cat("Writing to ", outfile, " .\n")
    write.table( param.df, file=outfile, row.names=FALSE, quote=FALSE )
}
EOF

for SIMDIR in ${BASEDIR}/test_*
do
    echo "Doing ${SIMDIR}."
    mkdir -p ${SIMDIR}/256x
    HTFILE=${SIMDIR}/256x/six-raster-list-hitting-times-full.tsv
    Rscript make-resistance-distances.R ../geolayers/multigrid/256x/crm_ 256x six-raster-list ${SIMDIR}/six-params.tsv analytic ${HTFILE}
    TRUTH_START=${BASEDIR}/inferred-six-params-true-start.tsv
    Rscript fit-logistic-model.R ../geolayers/multigrid/256x/crm_ 256x six-raster-list ${HTFILE} ${SIMDIR}/six-params.tsv ${TRUTH_START} &
    RANDOM_START=${BASEDIR}/inferred-six-params-random-start.tsv
    Rscript fit-logistic-model.R ../geolayers/multigrid/256x/crm_ 256x six-raster-list ${HTFILE} ${BASEPARAMS} ${RANDOM_START} &
    wait
done
