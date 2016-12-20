To compute depth, with output to files prefixed by `bamlist`:
```
ls *.bam > bamlist.txt
angsd -bam bamlist.txt -out 272torts_1e6_minmapq20minq30_map2kenro -uniqueOnly 1 -only_proper_pairs 1 -remove_bads 1 -doCounts 1 -nThreads 16 -minMapQ 20 -minQ 30
```
