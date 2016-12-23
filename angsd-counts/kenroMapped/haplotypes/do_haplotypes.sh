#!/bin/bash

LONG_SCAFS="
scaffold_288
scaffold_2062
scaffold_2603
scaffold_490
scaffold_2893
scaffold_2070
scaffold_54
scaffold_1314
scaffold_919
scaffold_1949
"

NS_SCAFS="
scaffold_272
scaffold_695
scaffold_1444
scaffold1752
scaffold3437
scaffold6830
scaffold8919
scaffold10886
scaffold18938
scaffold19403
scaffold19883
scaffold28862
scaffold38396
"

for x in $LONG_SCAFS $NS_SCAFS
do
    ../../../count-utils/plot-genotypes.R ../272torts_snp1e6_minmapq20minq30_map2kenro.counts.bin ../272torts_snp1e6_minmapq20minq30_map2kenro.pos.gz $x
done
