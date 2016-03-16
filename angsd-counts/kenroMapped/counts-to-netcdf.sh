#!/bin/bash

usage="
Convert a .counts.gz and .pos.gz file into a netcdf4 file.

Usage:
    $0 (counts file) (pos file) (output .nc file)
"

if [ $# -lt 3]
then
    echo $USAGE
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

NCDF_HEADER="
netcdf counts {

dimensions:
    indiv = 3, 
    base = 4,
    pos = unlimited;

variables:
    int    ${INDIVNAME}(indiv), ${POSNAME}(pos);
    string ${CHRNAME}(pos), ${BASESNAME}(base);
    int    ${COUNTNAME}(pos,indiv,base);  // last dim varies fastest

    ${INDIVNAME}:units = 'indiv';
    ${BASESNAME}:units = 'base';
    ${POSNAME}:units = 'snp';
    ${COUNTNAME}:units = 'reads';

data:
"
NCDF_FOOTER="
}
"

INDIVS=$(zcat $COUNTSFILE | head -n 1 | sed -e 's/_[ACGT]//g' | uniq | tr '\t' ',')
BASES=$(zcat $COUNTSFILE | head -n 1 | cut -f 1-4 | sed -e 's/[^\t_]*_//g' | tr '\t' ',')

echo "Writing to $OUTFILE "

( echo "$NCDF_HEADER";
  echo "${INDIVNAME} = ";
  echo "${INDIVS};";
  echo "${BASESNAME} = ";
  echo "${BASES};";
  echo "${CHRNAME} = ";
  zcat ${POSFILE} | tail -n +2 | cut -f 1 | tr '\n' ',' | sed -e 's/,$/;\n/';
  echo "${POSNAME} = ";
  zcat ${POSFILE} | tail -n +2 | cut -f 2 | tr '\n' ',' | sed -e 's/,$/;\n/';
  echo "${COUNTNAME} = ";
  zcat ${COUNTSFILE} | tail -n +2 | tr '\t' ',' | tr '\n' ',' | sed -e 's/,$/;\n/';
  echo ${NCDF_FOOTER}
  ) | ncgen -4 -o $OUTFILE 
