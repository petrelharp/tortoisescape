#!env Rscript

usage <- "Usage:
    Rscript get-coverage-stats.R (name of counts file) > stats.tsv
The counts file can be text, gzipped (.counts.gz) or binary (.counts.bin.gz, as output by 
counts-to-bin.R, one byte per int), in four-column _A _C _G _T format.

Outputs a tab-separated table of:
  nInd
  nAlleles
  total_coverage
  max_coverage
  min_coverage
  sd_coverage
  heterozygosity
to stdout.

A heterozygosity of -1 means that no sample had more than 1 read, so that heterozygosity 
could not be computed.
"

arglist <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(arglist)!=1) { stop(usage) }
countfile <- arglist[1]

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
    count.header <- scan(paste0(countfile,".header"),what="char")
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
attr(count.con,"nindivs") <- nindivs

# this returns an 3 x nsites matrix,
# then the rows are, for each allele:
#   nInd
#   nAlleles
#   total_A
#   total_C
#   total_G
#   total_T
#   max_coverage
#   min_coverage
#   sd_coverage
#   major_freq
#   heterozygosity
# columns of counts (transposed below) are of the form:
#   ind0_A  ind0_C  ind0_G  ind0_T  ind1_A  ind1_C  ind1_G  ind1_T  ...
# note that the output is transposed to what we might like,
# so that when we write out as a vector in chunks it will end up in the right order.
header.names <- c(
    "nInd",          # number of individuals with any reads
    "max_coverage",  # max total coverage across individuals
    "min_coverage",
    "sd_coverage",
    "total_A",       # number of A's
    "total_C",
    "total_G",
    "total_T",
    "heterozygosity" # mean heterozygosity of a random individual with coverage at least 2
    )
do_stats <- function (counts) {
    # assumes that the counts matrix has been "transposed" from the file,
    # i.e., read in with one row of the file equal to one column of counts
    # so that counts[4*(j-1)+k,] correponds to counts of the k'th nucleotide in the j-th individual

    # total coverage by individual
    coverage <- ( counts[1+4*(0:(nindivs-1)),] + counts[2+4*(0:(nindivs-1)),] 
                      + counts[3+4*(0:(nindivs-1)),] + counts[4+4*(0:(nindivs-1)),] )
    nInd <- colSums(coverage>0)
    max.coverage <- matrixStats::colMaxs(coverage)
    min.coverage <- matrixStats::colMins(coverage)
    sd.coverage <- matrixStats::colSds(coverage)
    # total coverage by allele
    totals <- sapply( 1:4, function (k) {
            colSums( counts[k+4*(0:(nindivs-1)),] )
        } )
    colnames(totals) <- paste0("total_",c("A","C","G","T"))
    # heterozygosity: for a single individual this is 1-sum_b n[b]*(n[b]-1)/n(n-1)
    homozygosity <- ( ( counts[1+4*(0:(nindivs-1)),] * (counts[1+4*(0:(nindivs-1)),]-1) 
                        + counts[2+4*(0:(nindivs-1)),] * (counts[2+4*(0:(nindivs-1)),]-1) 
                        + counts[3+4*(0:(nindivs-1)),] * (counts[3+4*(0:(nindivs-1)),]-1) 
                        + counts[4+4*(0:(nindivs-1)),] * (counts[4+4*(0:(nindivs-1)),]-1) )
                    / ( coverage*(coverage-1) ) )
    heterozygosity <- 1 - colMeans(homozygosity,na.rm=TRUE)
    heterozygosity[!is.finite(heterozygosity)] <- -1
    #
    return( rbind(
              nInd,
              max.coverage,
              min.coverage,
              sd.coverage,
              t(totals),
              heterozygosity
          ) )
}

# loop through the file
writeLines(paste0(header.names,collapse='\t'), headerfile)

nlines <- 0
while (length(counts <- read_fun(blocksize))>0) {
    dim(counts) <- c(4*nindivs,length(counts)/(4*nindivs))
    # rownames(counts) <- count.header
    nlines <- nlines+ncol(counts)
    cat("row ", nlines, "\n",file=stderr())
    stats <- do_stats(counts)
    write.table( t(stats), file=outcon, append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE )
    # for testing purposes:
    # if (nlines > 1e6) { break }
}
try(close(outcon))
try(close(count.con))


