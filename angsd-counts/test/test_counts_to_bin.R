## test counts-to-bin.R

countsfile <- "test.counts.gz"
posfile <- "test.pos.gz"
counts <- as.matrix(read.table(countsfile,header=TRUE))
pos <- read.table(posfile,header=TRUE)

binfile <- "test.bin"
nbytes <- 1
system(paste("Rscript ../counts-to-bin.R",countsfile,binfile,nbytes))

nindivs <- ncol(counts)/4

# to read lines from the result, NA'ing out counts above the maximum:
bincount <- file(binfile,open="rb")
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

if ( !all( c(
        all.equal( as.vector(counts), as.vector(read_bincounts(bincount,1:nrow(counts))) ),
        all.equal( as.vector(counts[3,]), as.vector(read_bincounts(bincount,3)) ),
        all.equal( as.vector(counts[7,]), as.vector(read_bincounts(bincount,7)) ),
        all.equal( as.vector(counts[3,]), as.vector(read_bincounts(bincount,3)) ) ) ) ) {
    stop("error!")
} else { cat("all good!\n") }

