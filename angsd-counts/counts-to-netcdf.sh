#!/bin/bash

USAGE="
Convert a .counts.gz and .pos.gz file into a netcdf4 file.

Usage:
    $0 (counts file) (pos file) (output .nc file)
"

if [ $# -lt 3 ]
then
    echo "$USAGE"
    exit 1
fi

COUNTSFILE=$1
POSFILE=$2
OUTFILE=$3

INDIVNAME="indiv"
POSNAME="position"
CHRNAME="chr"
BASESNAME="base"
COUNTNAME="count"

INDIVS=$(zcat $COUNTSFILE | head -n 1 | sed -e 's/_[ACGT]//g' | tr '\t' '\n' | uniq | sed -e 's/^/"/' -e 's/$/"/' | tr '\n' ',' | sed -e 's/,$//')
BASES=$(zcat $COUNTSFILE | head -n 1 | cut -f 1-4 | sed -e 's/[^\t_]*_//g' | sed -e 's/\t/", "/g' -e 's/^/"/' -e 's/$/"/')

NUM_INDIVS=$(echo $INDIVS |awk --field-separator="," "{ printf NF; } ")
NUM_BASES=$(echo $BASES | wc -w)  # paranoid a little?

NCDF_HEADER="
netcdf counts {

dimensions:
    indiv = ${NUM_INDIVS}, 
    base = ${NUM_BASES},
    pos = unlimited;

variables:
    int ${POSNAME}(pos);
    string ${CHRNAME}(pos), ${BASESNAME}(base), ${INDIVNAME}(indiv);
    int    ${COUNTNAME}(pos,indiv,base);  // last dim varies fastest

    ${INDIVNAME}:units = \"indiv\";
    ${BASESNAME}:units = \"base\";
    ${POSNAME}:units = \"snp\";
    ${COUNTNAME}:units = \"reads\";

data:
"
NCDF_FOOTER="
}
"


echo "Writing to $OUTFILE "

( echo "$NCDF_HEADER";
  echo "${INDIVNAME} = ";  # individuals
  echo "${INDIVS};";
  echo "${BASESNAME} = ";  # bases
  echo "${BASES};";
  echo "${CHRNAME} = ";  # chromosome
  zcat ${POSFILE} | tail -n +2 | cut -f 1 | sed -e 's/^/"/' -e 's/$/"/' | tr '\n' ',' | sed -e 's/,$//';
  echo ";"
  echo "${POSNAME} = ";  # position
  zcat ${POSFILE} | tail -n +2 | cut -f 2 | tr '\n' ',' | sed -e 's/,$//';
  echo ";"
  echo "${COUNTNAME} = ";  # counts
  zcat ${COUNTSFILE} | tail -n +2 | tr '\t' ',' | sed -e 's/$/,/g' | awk '{printf "%s",p} {p=$0 ORS} END{sub(/,.$/,"",p); print p}';
  echo ";"
  echo ${NCDF_FOOTER}
  ) | ncgen -4 -o $OUTFILE 


