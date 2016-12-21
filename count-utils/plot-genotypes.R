#!/usr/bin/Rscript

usage <- "
Usage:
   plot-genotypes.R (bincount file) (posfile) (scaffold name) [start:end]
produces a matrix plot of genotypes, with blue for alleles common in the north,
purple for alleles common in the south, red for heterozygotes, etc.

Example:
  plot-genotypes.R 272torts_snp1e6_minmapq20minq30_map2kenro.counts.bin 272torts_snp1e6_minmapq20minq30_map2kenro.pos.gz scaffold_272 
"

arglist <- commandArgs(TRUE)
if (length(arglist)<3) { stop(usage) }
bincountfile <- arglist[1]
posfile <- arglist[2]
scaffold <- arglist[3]
positions <- if (length(arglist)>3) { strsplit(arglist[4],":")[[1]] } else { c(0,Inf) }

# SNP info, and which lines the relevant SNPs are in
pos <- read.table(pipe(sprintf("zcat %s | grep -n '%s\\>' | sed -e 's/:/\t/'", posfile, scaffold)),
        header=FALSE,stringsAsFactors=FALSE)
names(pos) <- c("line","chr","pos","totDepth")
stopifnot(all(pos$chr==scaffold))
pos <- subset(pos,pos>=positions[1] & pos<=positions[2])

# sample info
tortdir <- gsub("tortoisescape.*","tortoisescape",getwd())
coord.obj <- load(file.path(tortdir,"tort_272_info","geog_coords.RData"))
coords <- get(coord.obj)
tort.ids <- row.names(coords)
nindivs <- length(tort.ids)
pcs <- read.csv(file.path(tortdir,"tort_272_info","pcs.csv"),header=TRUE,stringsAsFactors=FALSE)

# the counts themselves
bincount <- file(bincountfile,open="rb")
attr(bincount,"nindivs") <- nindivs
attr(bincount,"nbytes") <- 1
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

# read in the matrix of raw counts:
# this is (number of sites) x (4 * number of samples)
genotypes <- read_bincounts(bincount,pos$line)
