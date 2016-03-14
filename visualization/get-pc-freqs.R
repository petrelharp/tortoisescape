#!env Rscript

## Compute the inner product of PC1 with the counts of every site, for every base.
## Writes out a (nsites x 4) matrix of these inner products.

blocksize <- 1e5 # number of sites to read in at a time
do.text <- FALSE # write out in text? (if not, binary)
pc.num <- 2

# read in the site information
tortdir <- gsub("tortoisescape.*","tortoisescape",getwd())

# datadir <- file.path(dirname(tortdir),"angsd-counts")
# countfile <- file.path(datadir,"test.counts.gz")

# datadir <- file.path(tortdir,"angsd-counts")
# countfile <- file.path(datadir,"272torts_snp1e6_minmapq20minq30.counts.gz")
# maf <- read.table(file.path(datadir,"272torts_snp1e6_minmapq20minq30.mafs.gz"),header=TRUE)
# pos <- read.table(file.path(datadir,"272torts_snp1e6_minmapq20minq30.pos.gz"),header=TRUE)

datadir <- file.path(tortdir,"angsd-counts/kenroMapped")
countfile <- file.path(datadir,"272torts_snp1e6_minmapq20minq30_map2kenro.counts.gz")
# maf <- read.table(file.path(datadir,"272torts_snp1e6_minmapq20minq30_map2kenro.mafs.gz"),header=TRUE)
# pos <- read.table(file.path(datadir,"272torts_snp1e6_minmapq20minq30_map2kenro.pos.gz"),header=TRUE)

outsuffix <- if (do.text) { ".pc1counts.txt" } else { ".pc1counts.4bin" }
outfile <- gsub(".counts.gz",outsuffix,countfile)  # .4bin means binary, four columns

# the count file
count.con <- pipe(paste("zcat",countfile),open="r")
count.header <- scan(count.con,nlines=1,what="char")

###
# find pc-weighted counts by base
count.ids <- do.call(rbind,strsplit(count.header,"_"))
colnames(count.ids) <- c("angsd.id","base")
pcs <- read.csv(file.path(tortdir,"tort_272_info","pcs.csv"),header=TRUE,stringsAsFactors=FALSE)
# yes, they're in the same order (evan, 2/29/16)
pcs$angsd.id <- count.ids[count.ids[,2]=="A",1]

pcvec <- pcs[[paste0("PC",pc.num)]]

# this returns an 4 x nsites matrix,
# with the rows giving the inner product of the counts for each allele with PC1
# note that this is transposed to what we might like,
# so that we we write out as a vector in chunks it will end up in the right order.
nindivs <- nrow(count.ids)/4
do_pc_counts <- function (counts) {
    # assumes that the counts matrix has been "transposed" from the file,
    # i.e., read in with one row of the file equal to one column of counts
    # so that counts[4*(j-1)+k,] correponds to counts of the k'th nucleotide in the j-th individual
    do.call(rbind, lapply( 1:4, function (k) {
            as.vector(pcvec %*% counts[k+4*(0:(nindivs-1)),])
        } ) )
    # rownames(pc.counts) <- c("A","C","G","T")
}

# loop through the file
if (do.text) {
    writeLines("A\tC\tG\tT", outfile)
} else {
    outcon <- file(outfile, open="wb", raw=TRUE)
}
nlines <- 0
while (length(counts <- scan(count.con,nlines=blocksize))>0) {
    dim(counts) <- c(4*nindivs,length(counts)/(4*nindivs))
    # rownames(counts) <- count.header
    nlines <- nlines+ncol(counts)
    cat("row ", nlines, "\n")
    pc.counts <- do_pc_counts(counts)
    if (do.text) {
        write.table( t(pc.counts), file=outfile,  append=TRUE, sep='\t', row.names=FALSE, col.names=FALSE  )
    } else {
        writeBin( as.vector(pc.counts), con=outcon)
    }
    # for testing purposes:
    #   if (nlines > 1e6) { break }
}
try(close(outcon))
try(close(count.con))


if (FALSE) {  # to read it in
    pc.con <- pipe(paste("cat",outfile), open="rb")
    # do this multiple times to read in chunks
    pc.counts <- readBin(pc.con,what="numeric",n=4*blocksize)
    dim(pc.counts) <- c(4,length(pc.counts)/4)
    pc.counts <- t(pc.counts)
    colnames(pc.counts) <- c("A","C","G","T")
    # and then close
    close(pc.con)
}
