To compute depth, with output to files prefixed by `bamlist`:
```
ls ../*.bam | sort -V | grep -v "298\|299\|300\|301" > 272_bamlist.txt
angsd -bam 272_bamlist.txt -out 272torts_1e6_minmapq20minq30_map2kenro_hist -uniqueOnly 1 -only_proper_pairs 1 -remove_bads 1 -doCounts 1 -nThreads 16 -minMapQ 20 -minQ 30
```
