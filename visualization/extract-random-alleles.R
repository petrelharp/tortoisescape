#!env Rscript

nsamps <- 100 # number of alleles per frequency bin

# read in the site information
tortdir <- gsub("tortoisescape.*","tortoisescape",getwd())
datadir <- file.path(dirname(tortdir),"angsd-counts")
countfile <- file.path(datadir,"272torts_snp1e6_minmapq20minq30.counts.gz")

maf <- read.table(file.path(datadir,"272torts_snp1e6_minmapq20minq30.mafs.gz"),header=TRUE)
pos <- read.table(file.path(datadir,"272torts_snp1e6_minmapq20minq30.pos.gz"),header=TRUE)

# scaffolds
minscaf <- 1e3  # minimum number of snps per scaffold
scaf.counts <- table( pos$chr )
long.scafs <- names(scaf.counts)[scaf.counts>=minscaf]

# apply filters
mindepth <- 150 # minimum totalDepth
maxdepth <- 350 # maximum totalDepth
goodones <- ( ( pos$totDepth >= mindepth ) & ( pos$totDepth <= maxdepth ) )
# minind <- 250   # minimum nInd -- RESULTS IN NO SNPs
# goodones <- goodones & ( maf$nInd >= minind )
goodones <- goodones & ( pos$chr %in% long.scafs )

# select random sites
freq.breaks <- c(0,.05,(1:5)/10,1)  # frequency bins

site.list <- tapply( 
                      which(goodones),
                      cut(maf$knownEM[goodones],freq.breaks),
                      function (kk) {
                          sample(x=kk,size=min(length(kk),nsamps),replace=FALSE) }
                  )
site.info <- cbind( site=unlist(site.list), freqbin=rep(names(site.list),sapply(site.list,length)), pos[unlist(site.list),], maf[unlist(site.list),-(1:2)] )
site.info <- site.info[order(site.info$site),]

# extract those from the counts file
count.con <- gzfile(countfile,open="r")
count.header <- scan(count.con,nlines=1,what="char")
counts <- sapply( diff(c(1,1+site.info$site))-1, function (dk) { scan(count.con,skip=dk,nlines=1) } )
close(count.con)
rownames(counts) <- count.header

# write out to subdirectories
dirnames <- paste( "freq_", gsub("[^.0-9-]","",gsub(",","-",levels(site.info$freqbin))), sep='' )
for (k in 1:nlevels(site.info$freqbin)) {
    dir.create(dirnames[k])
    outbase <- file.path(dirnames[k],"random_sites")
    dothese <- (as.numeric(site.info$freqbin) == k)
    out.info.file <- paste(outbase,".info",sep='')
    out.count.file <- paste(outbase,".counts",sep='')
    write.table( site.info[dothese,], file=out.info.file, append=file.exists(out.info.file), sep='\t', row.names=FALSE, col.names=!file.exists(out.info.file) )
    write.table( t(counts[,dothese]), file=out.count.file,  append=file.exists(out.count.file), sep='\t', row.names=FALSE, col.names=!file.exists(out.count.file)  )
}

