#!/bin/bash
#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:sl230s:ppn=16
#PBS -l walltime=12:00:00
#PBS -l mem=24gb
#PBS -l vmem=120gb
#PBS -l pmem=1500mb

if [ ! -z ${PBS_O_WORKDIR-} ]
then
    source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh
    cd $PBS_O_WORKDIR
else
    DIRNAME=$1
    ALTNAME=$2
fi

if [ -z ${DIRNAME-} ] || [ -z ${ALTNAME-} ]
then
    echo "Usage:\
        evaluate-alternative.sh (name of directory) (name of alternative)
or
        qsub -vDIRNAME=(name of directory),ALTNAME=(name of alternative) evaluate-alternative.sh
"
    exit
fi

RMD=$(readlink -f reports/evaluate-alternative.Rmd)
OUTFILE="$DIRNAME/evaluate-${ALTNAME}.html"

# note that despite the setwd() below, the .Rmd does *not* behave as if in $DIRNAME.
ln -s -f $RMD $DIRNAME
R -e "setwd(\"${DIRNAME}\");library(knitr);alt.file=\"${ALTNAME}\";knit2html(\"evaluate-alternative.Rmd\",output=\"${OUTFILE}\");"

exit
