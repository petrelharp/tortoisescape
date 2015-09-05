#!/bin/bash

usage="Usage:
    $0 (maf file) (counts file) (min frequency) (max frequency) [number of sites]

Chooses random sites with MAF between [min] and [max]; outputs info and corresponding counts to stdout.

"

if [ $# -lt 4 ]
then
    echo $usage
    exit 1
fi

MAFFILE="$1"
COUNTFILE="$2"
MINFREQ="$3"
MAXFREQ="$4"
NSITES=${5:-1}  # defaults to 1

# 6th column is MAF
MAFINFO=$(zcat $MAFFILE | nl | awk -v minfreq=$MINFREQ -v maxfreq=$MAXFREQ '$6>minfreq && $6<maxfreq' | shuf -n $NSITES | sort)

SITES=$(echo "$MAFINFO" | cut -f 1)  # line numbers added by nl
AWKPAT=$(echo $(for x in $SITES; do echo "NR==$x ||"; done) | sed -e 's/||$//'; echo ";")

paste <(echo $MAFINFO | cut -f 2-) <(zcat $COUNTFILE | awk -f <(echo $AWKPAT))
