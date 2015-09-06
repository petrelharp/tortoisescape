#!/bin/bash

usage="Usage:
    $0  (min frequency) (max frequency) [number of sites] [scaffold]

Chooses random sites with MAF between [min] and [max]; 
optionally restricting to 'scaffold' (a regex) and at least a certain number of individuals (nInd);
outputs info and corresponding counts to stdout.

"

if [ $# -lt 4 ]
then
    echo $usage
    exit 1
fi

MAFFILE="272torts_snp1e6_minmapq20minq30.mafs.gz"
POSFILE="272torts_snp1e6_minmapq20minq30.pos.gz"
COUNTFILE="272torts_snp1e6_minmapq20minq30.counts.gz"

## all the info is:
# > paste <(zcat 272torts_snp1e6_minmapq20minq30.pos.gz) <(zcat 272torts_snp1e6_minmapq20minq30.mafs.gz | cut -f 3-) | nl
#  1  chr         pos  totDepth  major  minor  knownEM   pK-EM         nInd
#  2  scaffold_0  191  195       A      G      0.336019  0.000000e+00  132
#  3  scaffold_0  471  9         C      T      0.237115  1.676451e-07  8
#  4  scaffold_0  494  348       T      C      0.019132  4.440892e-16  204
#  5  scaffold_0  518  331       A      G      0.032015  0.000000e+00  193

MINFREQ="$1"
MAXFREQ="$2"
NSITES=${5:-1}  # defaults to 1
SCAFFOLD="${6:-.}"  # defaults to '.'
NIND="${7:-150}"  # defaults to 150

# 7th column is MAF (after nl)
MAFPAT="\$7>$MINFREQ && \$7<$MAXFREQ"
# 2nd is scaffold
if [ $# -ge 6 ]
then
    MAFPAT="$MAFPAT && \$2 ~ $SCAFFOLD"
fi

MAFINFO=$(paste <(zcat 272torts_snp1e6_minmapq20minq30.pos.gz) <(zcat 272torts_snp1e6_minmapq20minq30.mafs.gz | cut -f 3-) | nl | awk -f <(echo "$MAFPAT") | shuf -n $NSITES | sort)

SITES=$(echo "$MAFINFO" | cut -f 1)  # extract line numbers added by nl (note both files have a header)
AWKPAT=$(echo $(for x in $SITES; do echo "NR==$x ||"; done) | sed -e 's/||$/;/')

if [ -z "$SITES" ]
then
    echo "No sites found."
    exit 1
fi

paste <(echo "$MAFINFO" | cut -f 2-) <(zcat $COUNTFILE | awk -f <(echo "$AWKPAT"))
