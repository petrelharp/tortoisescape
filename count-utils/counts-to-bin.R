#!/usr/bin/env Rscript

usage <- "Convert read counts dumped by ANGSD to binary format.
Stores these in binary format, [nbytes] bytes each. The default is 1 byte (max value 255),
with values truncated at the maximum. Reads in roughly (chunksize) values at once (default 1e8).
Usage:\
    counts-to-bin.R (counts file) (output file) [chunksize] [nbytes]\
"

arglist <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if ( length(arglist) < 2 ) { stop(usage) }
countsfile <- arglist[1]
outfile <- arglist[2]
chunksize <- if (length(arglist)>2) { as.numeric(arglist[3]) } else { 1e8 }
nbytes <- if (length(arglist)>3) { as.integer(arglist[4]) } else { 1L }
maxcounts <- as.integer(256^nbytes-1)

if (file.exists(outfile)) {
    stop(sprintf("Output file %s already exists.", outfile))
}

#' Open for reading or writing
#'
#' Use to open stdin/stdout or process substitution things correctly
#'   from  http://stackoverflow.com/questions/15784373/process-substitution
#'
#' @param arg A text string: one of "-", "[/dev/]stdin", "[/dev/]stdout", or a file name (including a file descriptor, e.g. "/dev/fd3").
#'
#' @return A connection, from one of `stdout()`, `stdin()`, `fifo()`, or `file()`.
#'
#' @export
openread <- function(arg, open='r') {
    if (arg %in% c("-", "/dev/stdin","stdin")) {
       stdin()
    } else if (grepl("^/dev/fd/", arg)) {
       fifo(arg, open = open)
    } else {
       file(arg, open = open)
    }
}
openwrite <- function(arg, open='w') {
    if (arg %in% c("-", "/dev/stdout","stdout")) {
       stdout()
    } else if (grepl("^/dev/fd/", arg)) {
       fifo(arg, open = open)
    } else {
       file(arg, open = open)
    }
}

count.con <- openread(countsfile)
count.header <- scan(count.con,nlines=1,what="char")
count.ids <- do.call(rbind,strsplit(count.header,"_"))
colnames(count.ids) <- c("angsd.id","base")
nindivs <- nrow(count.ids)/4

headerfile <- paste0(outfile,".header")
cat( paste(count.header,collapse="\t"), "\n", sep='', file=headerfile )

# number of rows to read at once
chunklines <- floor(chunksize/length(count.header))

cat(paste("Writing to ", outfile, "\n"))
outcon <- openwrite(outfile,open="wb")

while (length(counts <- scan(count.con,what=integer(),nlines=chunklines))>0) {
    writeBin( pmin(counts,maxcounts), outcon, size=nbytes )
}

close(outcon)
try(close(count.con))

cat("All done writing to ", outfile, "\n")


if (FALSE) {
    # to read lines from the result, NA'ing out counts above the maximum:
    bincount <- file(outfile,open="rb")
    attr(bincount,"nindivs") <- nindivs
    attr(bincount,"nbytes") <- nbytes
    read_bincounts <- function (bincount,lines) {
        line.coords <- 4*(lines-1)*attr(bincount,"nbytes")*attr(bincount,"nindivs")
        output <- matrix( integer(4*length(lines)*attr(bincount,"nindivs")), nrow=length(lines) )
        for (k in seq_along(lines)) {
            seek(bincount,line.coords[k])
            output[k,] <- readBin( bincount, what=integer(), 
                                  n=4*attr(bincount,"nindivs"),
                                  size=attr(bincount,"nbytes"), 
                                  signed=(attr(bincount,"nbytes")>2) )
        }
        output[output>=256^attr(bincount,"nbytes")-1] <- NA
        return(output)
    }
}
