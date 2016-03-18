#!env Rscript

usage <- "Usage:\
    Rscript get-pc-counts.R (name of counts file) (pc number)\
\
The counts file can be text, gzipped (.counts.gz) or binary (.counts.bin, as output by counts-to-bin.R, one byte per int).\
Details:\
To compute the correlation and covariance of each site with a vector, say, PC1:\
If at a given site:\
    c(i) = major allele count in indiv i\
    n(i) = sample size, i.e. coverage of indiv i\
    p(i) = c(i)/n(i)\
    v(i) = vector of weights\
then we need these things:\
    A = sum_i p(i)\
    B = sum_i p(i) * v(i)\
    C = sum_{i : n(i)>0} v(i)\
    D = sum_{i : n(i)>0} v(i)^2\
    N = #{ i : n(i) > 0 }\
so we will estimate the covariance by\
    ( B - A * C / N ) / N\
and the correlation by\
    ( cov ) / sqrt( ( ( A - A^2 / N ) / N ) * ( ( D - C^2 / N ) / N ) )\
        = ( B * N - A * C ) / sqrt( ( A * N - A^2 ) * ( D * N - C^2 ) )\
"

arglist <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(arglist)!=2) { stop(usage) }
countfile <- arglist[1]
pc.num <- as.integer(arglist[2])

# pc.num <- 1

# datadir <- file.path(tortdir,"angsd-counts/kenroMapped")
# countfile <- file.path(datadir,"272torts_snp1e6_minmapq20minq30_map2kenro.counts.gz")
# maf <- read.table(file.path(datadir,"272torts_snp1e6_minmapq20minq30_map2kenro.mafs.gz"),header=TRUE)
# pos <- read.table(file.path(datadir,"272torts_snp1e6_minmapq20minq30_map2kenro.pos.gz"),header=TRUE)

blocksize <- 1e5 # number of sites to read in at a time
do.text <- FALSE # write out in text? (if not, binary)

# read in the site information
tortdir <- gsub("tortoisescape.*","tortoisescape",getwd())
outsuffix <- sprintf(if (do.text) { ".pc%dcounts.txt" } else { ".pc%dcounts.5bin" }, pc.num)
outfile <- gsub(".counts.gz",outsuffix,countfile)  # .5bin means binary, five columns
headerfile <- if ( do.text ) { outfile } else { paste0(outfile,".header") }

# the count file
if (grepl(".counts.gz$",countfile)) {
    count.con <- pipe(paste("zcat",countfile),open="r")
    count.header <- scan(count.con,nlines=1,what="char")
    read_fun <- function (blocksize) { scan(count.con,nlines=blocksize) }
} else if (grepl(".counts.bin$",countfile)) {
    count.con <- file(countfile,open="rb")
    count.header <- scan(paste0(countfile,".header"),what="char")
    attr(count.con,"nbytes") <- 1
    read_fun <- function (blocksize) {
        readBin( bincount, what=integer(),
                  n=4*attr(bincount,"nindivs")*blocksize,
                  size=attr(bincount,"nbytes"),
                  signed=(attr(bincount,"nbytes")>2) )
    }
} else { stop(paste("Counts file", countfile, "not a recognized format.")) }

## count info
count.ids <- do.call(rbind,strsplit(count.header,"_"))
colnames(count.ids) <- c("angsd.id","base")
nindivs <- nrow(count.ids)/4
attr(count.con,"nindivs") <- nindivs

###
# find pc-weighted counts by base
pcs <- read.csv(file.path(tortdir,"tort_272_info","pcs.csv"),header=TRUE,stringsAsFactors=FALSE)
# yes, they're in the same order (evan, 2/29/16)
pcs$angsd.id <- count.ids[count.ids[,2]=="A",1]

pcvec <- pcs[[paste0("PC",pc.num)]]

# this returns an 3 x nsites matrix,
# then the rows are, for each allele:
#   A = sum_i p(i)         -- sum of major allele frequencies
#   B = sum_i p(i) * v(i)  -- inner product of major allele frequencies with PC1
#   C = sum_{i : n(i)>0} v(i)
#   D = sum_{i : n(i)>0} v(i)^2
#   N = #{ i : n(i) > 0 }  -- number of sampled individuals
# columns of counts (transposed below) are of the form:
#   ind0_A  ind0_C  ind0_G  ind0_T  ind1_A  ind1_C  ind1_G  ind1_T  ...
# note that the output is transposed to what we might like,
# so that we we write out as a vector in chunks it will end up in the right order.
do_pc_counts <- function (counts) {
    # assumes that the counts matrix has been "transposed" from the file,
    # i.e., read in with one row of the file equal to one column of counts
    # so that counts[4*(j-1)+k,] correponds to counts of the k'th nucleotide in the j-th individual
    coverage <- ( counts[1+4*(0:(nindivs-1)),] + counts[2+4*(0:(nindivs-1)),] 
                      + counts[3+4*(0:(nindivs-1)),] + counts[4+4*(0:(nindivs-1)),] )
    totals <- sapply( 1:4, function (k) {
            colSums( counts[k+4*(0:(nindivs-1)),] )
        } )
    # names(totals) <- rownames(acgt.prod) <- c("A","C","G","T")
    max.counts <- pmax( totals[,1], totals[,2], totals[,3], totals[,4] )
    major.col <- ifelse( totals[,1]==max.counts, 1,
                        ifelse( totals[,2]==max.counts, 2,
                            ifelse( totals[,3]==max.counts, 3, 4 ) ) )
    # major allele frequencies
    major.freqs <- ( counts[ cbind( as.vector(outer(4*(0:(nindivs-1)),major.col,"+")), rep(1:ncol(counts),each=nindivs) ) ]
                    / as.vector(coverage) )
    major.freqs[!is.finite(major.freqs)] <- 0
    dim(major.freqs) <- c(nindivs,ncol(counts))
    return( rbind(
              colSums(major.freqs),                 # this is A
              as.vector(pcvec %*% major.freqs),     # this is B
              as.vector(pcvec %*% (coverage>0)),    # this is C
              as.vector(pcvec^2 %*% (coverage>0)),  # this is D
              colSums( coverage>0 )                 # this is N
          ) )
}

# loop through the file
writeLines("sum_freq\tfreq_prod\tsum_weights\tsum_weights_sq\tnum_nonzero\n", headerfile)
if (!do.text) {
    outcon <- file(outfile, open="wb", raw=TRUE)
}

nlines <- 0
while (length(counts <- read_fun(blocksize))>0) {
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
    pc.header <- scan(headerfile,what='char')
    pc.con <- pipe(paste("cat",outfile), open="rb")
    pc.header <-  scan(headerfile,what='char')
    # do this multiple times to read in chunks
    read_chunk <- function (blocksize) {
        # do this multiple times to read in chunks
        pc.counts <- readBin(pc.con,what="numeric",n=length(pc.header)*blocksize)
        dim(pc.counts) <- c(length(pc.header),length(pc.counts)/length(pc.header))
        pc.counts <- t(pc.counts)
        colnames(pc.counts) <- pc.header
        return(pc.counts)
    }
    # and then close
    close(pc.con)
}
