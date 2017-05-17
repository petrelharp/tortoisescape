#!/usr/bin/env Rscript

usage <- "Usage:\
    get-coverage-by-scaffold.R (name of counts file) (min total coverage) (max total coverage) (output file) [blocklines]\
\
Writes out a table of the form\
\
scaffold_name  start  end     num_sites  id1  id2  id3  ...\
scaffold_1     102    284272  1023       3    2    0    ...\
\
in which one row correponds to one scaffold, and records the positions of the first and last sites,\
the number of sites, and for each individual their total coverage on that scaffold.\
Only sites with total coverage between (min total coverage) and (max total coverage) are considered.\
Reads in (blocklines) lines at a time (default: 1e5)\
\
Assumes that if the counts file is\
  the-counts-file.counts.bin.gz
then there is also\
  the-counts-file.counts.header\
  the-counts-file.pos.gz\
containing information about columns and rows, respectively.
"

arglist <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }

if (length(arglist)<4) { stop(usage) }
countfile <- arglist[1]
min_coverage <- as.integer(arglist[2])
max_coverage <- as.integer(arglist[3])
outfile <- arglist[4]
# number of lines to read in at a time
blocklines <- if (length(arglist) > 4) { as.integer(arglist[5]) } else { 1e5 }

if (file.exists(outfile)) {
    stop(sprintf("File %s already exists.", outfile))
}

do.text <- TRUE # write out in text? (if not, binary)

posfile <- paste0(gsub("[.]counts[.].*","",countfile),".pos.gz")
if (!file.exists(posfile)) { stop(sprintf("Cannot find position file %s.", posfile)) }
# pos.con <- gzfile(posfile, open="r")
# pos.header <- scan(pos.con, nlines=1, what='')
# stopifnot(all(pos.header==c("chr", "pos", "totDepth"))) 
# pos_fun <- function (blocklines) { read.table(pos.con, nrows=blocklines) }
full.pos <- read.table(posfile, header=TRUE, stringsAsFactors=FALSE)
pos_fun <- function (blocklines) {
    out <- full.pos[pos.counter:min(nrow(full.pos),(pos.counter+blocklines-1)),] 
    assign("pos.counter", pos.counter+blocklines, parent.env(environment()))
    return(out)
}
environment(pos_fun) <- new.env()
assign("pos.counter", 1, environment(pos_fun))

# the count file
if (grepl(".counts.gz$",countfile)) {
    count.con <- gzfile(countfile,open="r")
    count.header <- scan(count.con,nlines=1,what="char")
    read_fun <- function (blocklines) { scan(count.con,nlines=blocklines) }
} else if (grepl("counts.bin",countfile)) {
    count.con <- if (!grepl(".counts.bin.gz$",countfile)) { file(countfile,open="rb") } else { gzfile(countfile,open="rb") }
    count.header <- scan(paste0(gsub(".gz$","",countfile),".header"),what="char")
    attr(count.con,"nbytes") <- 1
    read_fun <- function (blocklines) {
        readBin( count.con, what=integer(),
                  n=4*attr(count.con,"nindivs")*blocklines,
                  size=attr(count.con,"nbytes"),
                  signed=(attr(count.con,"nbytes")>2) )
    }
} else { stop(paste("Counts file", countfile, "not a recognized format.")) }

## count info
count.ids <- do.call(rbind,strsplit(count.header,"_"))
colnames(count.ids) <- c("angsd.id","base")
nindivs <- nrow(count.ids)/4
attr(count.con,"nindivs") <- nindivs
indiv.ids <- unique(count.ids[,1])
id_factor <- factor(count.ids[,"angsd.id"], levels=indiv.ids)

if (do.text) {
    headerfile <- outfile
} else {
    headerfile <- gsub(".gz$", ".header", outfile)
    outcon <- file(outfile, open="wb", raw=TRUE)
}

# loop through the file
header_names <- c(c("scaffold_name", "start", "end", "num_sites"), indiv.ids)
writeLines(paste(header_names, collapse="\t"), headerfile)


write_bit <- function (bit, nums) {
    out <- cbind(bit, nums)
    if (do.text) {
        write.table(out, file=outfile,  append=TRUE, sep='\t', 
                    row.names=FALSE, col.names=FALSE  )
    } else {
        writeBin( out, con=outcon)
    }
    # for testing purposes:
    #   if (nlines > 1e6) { break }
}

nlines <- 0
last_bit <- data.frame(scafs=character(0), starts=integer(0), ends=integer(0), num_sites=integer(0))
last_nums <- matrix(0, nrow=0, ncol=nindivs)
As <- 1+4*seq(0,length.out=nindivs)

while (length(counts <- read_fun(blocklines))>0) {
    dim(counts) <- c(4*nindivs,length(counts)/(4*nindivs))
    coverages <- counts[As,,drop=FALSE]
    for (k in 1:3) {
        coverages <- coverages + counts[As+k,,drop=FALSE]
    }
    total_coverages <- colSums(coverages)
    use_these <- (total_coverages < max_coverage) & (total_coverages > min_coverage)
    pos <- pos_fun(blocklines)
    if (!any(use_these)) { next }
    coverages <- coverages[,use_these,drop=FALSE]
    pos <- pos[use_these,,drop=FALSE]
    pos$chr <- factor(pos$chr, levels=unique(pos$chr))
    next_bit <- data.frame(
                    scafs=levels(pos$chr),
                    starts=tapply(pos$pos, pos$chr, min),
                    ends=tapply(pos$pos, pos$chr, max),
                    num_sites=tapply(pos$pos, pos$chr, length),
                    stringsAsFactors=FALSE)
    next_nums <- do.call(rbind, tapply(1:ncol(coverages), pos$chr, 
                                       function (kk) { rowSums(coverages[,kk,drop=FALSE]) },
                                   simplify=FALSE))
    if (nrow(last_bit) > 0) {
       if (next_bit$scafs[1] == last_bit$scafs) {
           next_bit$starts[1] <- last_bit$starts
           next_bit$num_sites[1] <- next_bit$num_sites[1] + last_bit$num_sites
           next_nums[1,] <- next_nums[1,] + last_nums
       } else {
           write_bit(last_bit, last_nums)
       }
    }
    if (nrow(next_bit)>1) {
        write_bit(next_bit[-nrow(next_bit),,drop=FALSE], next_nums[-nrow(next_bit),,drop=FALSE])
    }
    last_bit <- next_bit[nrow(next_bit),,drop=FALSE]
    last_nums <- next_nums[nrow(next_bit),,drop=FALSE]
}
write_bit(last_bit, last_nums)

if (!do.text) { try(close(outcon)) }
try(close(count.con))

