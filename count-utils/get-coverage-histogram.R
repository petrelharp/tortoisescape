#!/usr/bin/env Rscript

usage <- "Usage:
    Rscript get-coverage-histogram.R (name of counts file) (maximum count) > coverage.tsv
The counts file can be text, gzipped (.counts.gz) or binary (.counts.bin.gz, as output by 
counts-to-bin.R, one byte per int), in four-column _A _C _G _T format.

This outputs the per-sample coverage histograms to standard out.
"

arglist <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(arglist)!=2) { stop(usage) }
countfile <- arglist[1]
maxcount <- as.numeric(arglist[2])

# for colMax
library(matrixStats)

## output file
outcon <- headerfile <- stdout()

blocksize <- 1e5 # number of sites to read in at a time

# the count file
if (grepl(".counts.gz$",countfile)) {
    count.con <- gzfile(countfile,open="r")
    count.header <- scan(count.con,nlines=1,what="char")
    read_fun <- function (blocksize) { scan(count.con,nlines=blocksize) }
} else if (grepl(".counts.bin",countfile)) {
    count.con <- if (!grepl(".counts.bin.gz$",countfile)) { file(countfile,open="rb") } else { gzfile(countfile,open="rb") }
    count.header <- scan(paste0(gsub(".gz$","",countfile),".header"),what="char")
    attr(count.con,"nbytes") <- 1
    read_fun <- function (blocksize) {
        readBin( count.con, what=integer(),
                  n=4*attr(count.con,"nindivs")*blocksize,
                  size=attr(count.con,"nbytes"),
                  signed=(attr(count.con,"nbytes")>2) )
    }
} else { stop(paste("Counts file", countfile, "not a recognized format.")) }

## count info
count.ids <- do.call(rbind,strsplit(count.header,"_"))
colnames(count.ids) <- c("angsd.id","base")
nindivs <- nrow(count.ids)/4
indiv.ids <- unique(count.ids[,1])
attr(count.con,"nindivs") <- nindivs

header.names <- indiv.ids
do_stats <- function (counts) {
    # assumes that the counts matrix has been "transposed" from the file,
    # i.e., read in with one row of the file equal to one column of counts
    # so that counts[4*(j-1)+k,] correponds to counts of the k'th nucleotide in the j-th individual

    # total coverage by individual
    coverage <- ( counts[1+4*(0:(nindivs-1)),] + counts[2+4*(0:(nindivs-1)),] 
                      + counts[3+4*(0:(nindivs-1)),] + counts[4+4*(0:(nindivs-1)),] )
    coverage <- pmin(coverage,maxcount)
    dim(coverage) <- c(dim(counts)[1]/4,dim(counts)[2])
    d.hists <- apply(coverage,1,tabulate,nbins=maxcount)  # tabulate does *positive* integers
    # also return in "transposed" format
    return( rbind( ncol(coverage)-colSums(d.hists), d.hists) )
}

# loop through the file
writeLines(paste0(header.names,collapse='\t'), headerfile)

## the table of results
hists <- matrix(0,nrow=1+maxcount,ncol=nindivs)

nlines <- 0
while (length(counts <- read_fun(blocksize))>0) {
    dim(counts) <- c(4*nindivs,length(counts)/(4*nindivs))
    # rownames(counts) <- count.header
    nlines <- nlines+ncol(counts)
    cat("row ", nlines, "\n",file=stderr())
    hists <- hists + do_stats(counts)
    # for testing purposes:
    # if (nlines > 1e6) { break }
}

write.table( hists, file=outcon, append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE )

try(close(count.con))


