#!/bin/bash

MINQ=30
MINMAPQ=20
MAXDEPTH=100

USAGE="
Uses samtools to find the histogram of coverage, 
using a minimum quality score of $MINQ
and a minimum mapping quality of $MINMAPQ,
and truncated a maximum depth of $MAXDEPTH.
Usage:
    $0 (bamfile) > (coverage histogram)
"

if [ $# -lt 1 ]
then
    echo "$USAGE"
    exit 1
fi

# samtools depth -a output is like :
# (chrom)     pos depth
# scaffold_0	1	0
# scaffold_0	2	0
# scaffold_0	3	0
# scaffold_0	4	0
# scaffold_0	5	0
# scaffold_0	6	0

samtools depth -a -q $MINQ -Q $MINMAPQ $1 | awk 'BEGIN { for (i=0;i<='$MAXDEPTH';i++) n[i]=0 } { if ($3 > '$MAXDEPTH') { n['$MAXDEPTH']++; } else { n[$3]++; } } END {for (i=0;i<='$MAXDEPTH';i++) print i,n[i]}'

