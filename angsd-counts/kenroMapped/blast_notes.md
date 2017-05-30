
To find CIRBP:

Downloaded Xenopus CIRBP from NCBI Reference Sequence: NC_030677.1 to `CIRBP.fa`, which has 6181 bases.

Set up BLAST:
```
# make blast db to blast against
makeblastdb -in GopAga1.0.softmasked.fa -out GopAga1.0.softmasked.blastdb -input_type fasta -dbtype nucl
# blast
blastn -query CIRBP.fa -task blastn -db GopAga1.0.softmasked.blastdb -out CIRBP.blasted
blastn -query CIRBP.fa -task dc-megablast -db GopAga1.0.softmasked.blastdb -out CIRBP.dc-megablast
```

Results:

- `scaffold23396` - 90% identity, 
- `scaffold11893` - 83% identity
- `scaffold_510` - 83% identity but shorter


```r
scafs <- c("scaffold23396", "scaffold11893", "scaffold_510")

pos <- read.table(file.path(datadir,"272torts_snp1e6_minmapq20minq30_map2kenro.pos.gz"),header=TRUE,stringsAsFactors=FALSE)

for (scaf in scafs) {
    scafdir <- paste0(scaf,"_maps")
    dir.create(scafdir,showWarnings=FALSE)
    thisdir <- file.path( scafdir, "centered" )
    dir.create(thisdir,showWarnings=FALSE)
    scaf.lines <- which( pos$chr==scaf )
    do.lines <- scaf.lines[floor(length(scaf.lines)/2)+(1:100)]
    png( filename=file.path(thisdir,"snp_map%03d.png"), 
            width=5*144, height=5*144, pointsize=10, res=144 )
    map_alleles( do.lines )
    dev.off()
    system(paste("cd", thisdir, ";",
                 file.path(tortdir,"visualization/make-slideshow.sh"),
                 "*.png > index.html" ) )
}
```
