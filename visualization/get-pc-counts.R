#!env Rscript

## Compute the inner product of PC1 with the counts of every site, for every base.
## Writes out a (nsites x 4) matrix of these inner products.

blocksize <- 1e5 # number of sites to read in at a time
do.text <- FALSE # write out in text? (if not, binary)

# read in the site information
tortdir <- file.path(gsub("tortoises.*","tortoises",getwd()),"tortoisescape")

# datadir <- file.path(dirname(tortdir),"angsd-counts")
# countfile <- file.path(datadir,"272torts_snp1e6_minmapq20minq30.counts.gz")
# maf <- read.table(file.path(datadir,"272torts_snp1e6_minmapq20minq30.mafs.gz"),header=TRUE)
# pos <- read.table(file.path(datadir,"272torts_snp1e6_minmapq20minq30.pos.gz"),header=TRUE)

datadir <- file.path(dirname(tortdir),"angsd-counts/kenroMapped")
countfile <- file.path(datadir,"272torts_snp1e6_minmapq20minq30_map2kenro.counts.gz")
# maf <- read.table(file.path(datadir,"272torts_snp1e6_minmapq20minq30_map2kenro.mafs.gz"),header=TRUE)
# pos <- read.table(file.path(datadir,"272torts_snp1e6_minmapq20minq30_map2kenro.pos.gz"),header=TRUE)

outfile <- gsub(".counts.gz",".pccounts.4bin",countfile)  # .4bin means binary, four columns

# the count file
count.con <- gzfile(countfile,open="r")
count.header <- scan(count.con,nlines=1,what="char")

###
# find pc-weighted counts by base
count.ids <- do.call(rbind,strsplit(count.header,"_"))
colnames(count.ids) <- c("angsd.id","base")
pcs <- read.csv(file.path(tortdir,"tort_272_info","pcs.csv"),header=TRUE,stringsAsFactors=FALSE)
# yes, they're in the same order (evan, 2/29/16)
pcs$angsd.id <- count.ids[count.ids[,2]=="A",1]

pcvec <- pcs$PC1

# this returns an nsites x 4 matrix,
# with the columns giving the inner product of the counts for each allele
#   with PC1
nindivs <- nrow(count.ids)/4
do_pc_counts <- function (counts) {
    # assumes that the counts matrix has been "transposed" from the file,
    # i.e., read in with one row of the file equal to one column of counts
    sapply( 1:4, function (k) {
            counts[,k-1+(1:nindivs)] %*% pcvec
        } )
    # colnames(pc.counts) <- c("A","C","G","T")
}

# loop through the file
if (do.text) {
    writeLines("A\tC\tG\tT", outfile)
} else {
    outcon <- file(outfile, open="wb", raw=TRUE)
}
skip <- 1
while (length(counts <- scan(countfile,nmax=4*nindivs*blocksize,skip=skip))>0) {
    dim(counts) <- c(length(counts)/(4*nindivs),4*nindivs)
    # rownames(counts) <- count.header
    skip <- skip+nrow(counts)
    cat("row ", skip, "\n")
    pc.counts <- do_pc_counts(counts)
    if (do.text) {
        write.table( as.vector(pc.counts), file=outfile,  append=TRUE, sep='\t', row.names=FALSE, col.names=FALSE  )
    } else {
        writeBin(as.vector(pc.counts), con=outcon)
    }
    # for testing purposes:
    #   if (skip > 1e6) { break }
}
if (!do.text) { close(outcon) }


if (FALSE) {  # to read it in
    pc.con <- pipe(paste("cat",outfile), open="rb")
    # do this multiple times to read in chunks
    pc.counts <- readBin(pc.con,what="numeric",n=4*blocksize)
    dim(pc.counts) <- c(length(pc.counts)/4,4)
    # and then close
    close(pc.con)
}
