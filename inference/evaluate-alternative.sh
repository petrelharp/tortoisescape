#!/bin/bash
#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:sl230s:ppn=16
#PBS -l walltime=12:00:00
#PBS -l mem=24gb
#PBS -l vmem=64gb
#PBS -l pmem=1500mb

if [ ! -z ${PBS_O_WORKDIR-} ]
then
    source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh
    cd $PBS_O_WORKDIR
else
    DIRNAME=$1
    ALTNAME=$2
    SUMMARY=$3
fi

if [ -z ${DIRNAME-} ] || [ -z ${ALTNAME-} ] || [ -z ${SUMMARY-} ]
then
    echo "Usage:\
        evaluate-alternative.sh (name of directory) (name of alternative) (name of summary directory)
e.g.
        evaluate-alternative.sh nussear-south-flat/habitat-gt-zero alt_pref_pda habitat-only
or
        qsub -vDIRNAME=(name of directory),ALTNAME=(name of alternative),SUMMARY=(name of summary) evaluate-alternative.sh
"
    exit
fi

RMD=$(readlink -f reports/evaluate-alternative.Rmd)
OUTFILE="evaluate-${ALTNAME}.html"

# note that despite the setwd() below, the .Rmd does *not* behave as if in $DIRNAME.
ln -s -f $RMD $DIRNAME
R -e "setwd(\"${DIRNAME}\");library(knitr);alt.file=\"${ALTNAME}\";summary.dir=\"${SUMMARY}\";knit2html(\"evaluate-alternative.Rmd\",output=\"${OUTFILE}\");"

exit
