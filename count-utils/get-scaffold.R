#!/usr/bin/env Rscript

usage <- "Usage:\
    get-scaffold.R (name of counts file) (name of scaffold) > output\
Assumes that if the counts file is\
  the-counts-file.counts.bin.gz
then there is also\
  the-counts-file.counts.header\
  the-counts-file.pos.gz\
containing information about columns and rows, respectively.
"


arglist <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }

if (length(arglist)!=2) { stop(usage) }
countfile <- arglist[1]
scaffold <- arglist[2]

posfile <- paste0(gsub("[.]counts[.].*","",countfile),".pos.gz")
if (!file.exists(posfile)) { stop(sprintf("Cannot find position file %s", posfile)) }
awkscript <- sprintf("zcat %s | awk 'BEGIN {z=0} /%s/ {if (z==0) print FNR; z=1} !/%s/ {if (z==1) { print FNR; exit} } END { print FNR+1 }'", posfile, scaffold, scaffold)
# minus one for the header
scaf_lines <- scan(pipe(awkscript)) - 1

# the count file
if (grepl(".counts.gz$",countfile)) {
    count.con <- gzfile(countfile,open="r")
    count.header <- scan(count.con,nlines=1,what="char")
    read_fun <- function (start, end) { scan(count.con, skip=start-1, nlines=end-start) }
} else if (grepl("counts.bin",countfile)) {
    count.con <- if (!grepl(".counts.bin.gz$",countfile)) { file(countfile,open="rb") } else { gzfile(countfile,open="rb") }
    count.header <- scan(paste0(gsub(".gz$","",countfile),".header"),what="char")
    attr(count.con,"nbytes") <- 1
    # Read the lines from start to end, inclusive.
    read_fun <- function (start, end) {
        line_length <- 4*attr(count.con,"nindivs")
        # note that seek on gzfiles starts at zero and only moves forwards
        #  --> zero-based index
        seek(count.con, where=(start-1)*line_length, rw="r")
        readBin( count.con, what=integer(),
                  n=line_length*(end-start),
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

out <- read_fun(scaf_lines[1], scaf_lines[2])
# dim(out) <- c(line_length, end-start+1)

write(out, file=stdout(), ncolumns=4*nindivs)
