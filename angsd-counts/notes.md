To compute depth, with output to files prefixed by `bamlist`:

```
# get the list of bams
ls ../*.bam | sort -V | grep -v "298\|299\|300\|301" > 272_bamlist.txt

# and get a list of the scaffolds at least 10kb:
cat GopAga1.0.softmasked.scafflength | awk ' $2 >= 10000 { print; } ' > GopAga1.0.softmasked.scafflength.10kb
tail -n +2 GopAga1.0.softmasked.scafflength | awk ' $2 >= 10000 { print $1":"; } ' > scafnames_10kb.txt

split -l 120 scafnames_10kb.txt coverage/scafnames_10kb_sub_

# note that sequential-joblist_nproc_16.pbs changes directory first
# and needs more than 60 GB (or else segfault)

rm coverage/all_scafjobs.sh*

for sn in $(ls coverage/scafnames_10kb_sub_* | sed -e 's/coverage.//g')
do
    echo "angsd -bam 272_bamlist.txt -rf coverage/$sn -out coverage/272torts_1e6_minmapq20minq30_map2kenro_${sn}_hist -doDepth 1 -uniqueOnly 1 -only_proper_pairs 1 -remove_bads 1 -doCounts 1 -nThreads 16 -minMapQ 20 -minQ 30 -maxDepth 1200" >>coverage/all_scafjobs.sh
done

./parallel-sequential-joblist_nproc_16.sh 6 coverage/all_scafjobs.sh
```

And, combine the outputs:
```r
sample.files <- list.files("coverage","*.depthSample",full.names=TRUE)
out <- matrix(0, nrow=272, ncol=1201)
for (sf in sample.files) {
    x <- read.table(sf)
    stopifnot(all(dim(x)==dim(out)))
    out <- out + x
}
write.table(out, file="272torts_1e6_minmapq20minq30_map2kenro_scafnames_10kb_hist.depthSample",row.names=FALSE,col.names=FALSE)

global.files <- list.files("coverage","*.depthGlobal",full.names=TRUE)
out <- numeric(1201)
for (sf in global.files) {
    x <- read.table(sf)
    stopifnot(length(x)==length(out))
    out <- out + x
}
write.table(out, file="272torts_1e6_minmapq20minq30_map2kenro_scafnames_10kb_hist.depthGlobal",row.names=FALSE,col.names=FALSE)

```


Or, to do more of the contigs:

```
# this'd do everything but it takes FOREVER because the last (short) chroms take FOREVER for some reason
#   angsd -bam 272_bamlist.txt -out 272torts_1e6_minmapq20minq30_map2kenro_hist -doDepth 1 -uniqueOnly 1 -only_proper_pairs 1 -remove_bads 1 -doCounts 1 -nThreads 16 -minMapQ 20 -minQ 30
# so split it out into two jobs:
zcat 272torts_snp1e6_minmapq20minq30_map2kenro.pos.gz | cut -f 1 | grep scaffold | uniq | sed -e 's/$/:/' > scaffold_scaffold_names.txt
zcat 272torts_snp1e6_minmapq20minq30_map2kenro.pos.gz | cut -f 1 | tail -n +2 | grep -v scaffold | uniq | sed -e 's/$/:/' > scaffold_other_names.txt
# also it only ever uses like 2 CPU for some reason
echo "cd $PWD; angsd -bam 272_bamlist.txt -rf scaffold_scaffold_names.txt -out 272torts_1e6_minmapq20minq30_map2kenro_scaffoldscaffold_hist -doDepth 1 -uniqueOnly 1 -only_proper_pairs 1 -remove_bads 1 -doCounts 1 -nThreads 4 -minMapQ 20 -minQ 30" | qsub -q cmb -l nodes=1:ppn=4 -l mem=64gb -l walltime=120:00:00
echo "cd $PWD; angsd -bam 272_bamlist.txt -rf scaffold_other_names.txt -out 272torts_1e6_minmapq20minq30_map2kenro_scaffoldother_hist -doDepth 1 -uniqueOnly 1 -only_proper_pairs 1 -remove_bads 1 -doCounts 1 -nThreads 4 -minMapQ 20 -minQ 30" | qsub -q cmb -l nodes=1:ppn=4 -l mem=64gb -l walltime=120:00:00
```


