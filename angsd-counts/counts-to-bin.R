#!env Rscript

usage <- "Convert read counts dumped by ANGSD to binary format.
Stores these in binary format, [nbytes] bytes each. The default is 1 byte (max value 255),
with values truncated at the maximum.\
Usage:\
    counts-to-bin.R (counts file) (output file) [nbytes]\
"

arglist <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if ( length(arglist) < 2 ) { stop(usage) }
countsfile <- arglist[1]
outfile <- arglist[2]
nbytes <- if (length(arglist)>2) { as.integer(arglist[3]) } else { 1L }
maxcounts <- as.integer(256^nbytes-1)

count.con <- pipe(paste("zcat",countsfile),open="r")
count.header <- scan(count.con,nlines=1,what="char")
count.ids <- do.call(rbind,strsplit(count.header,"_"))
colnames(count.ids) <- c("angsd.id","base")
nindivs <- nrow(count.ids)/4

headerfile <- paste0(outfile,".header")
cat( paste(count.header,collapse="\t"), "\n", sep='', file=headerfile )

# number of values to read at once
chunksize <- 1e8
# number of rows to read at once
chunklines <- floor(chunksize/length(count.header))

cat(paste("Writing to ", outfile, "\n"))
outcon <- file(outfile,open="wb")

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
