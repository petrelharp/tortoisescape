#!/usr/bin/env Rscript
library(testthat)

set.seed(23)

## angsd output looks like this:
# ind0_A	ind0_C	ind0_G	ind0_T	ind1_A	ind1_C	ind1_G	ind1_T	ind2_A	ind2_C	ind2_G	ind2_T
# 0	0	0	0	0	0	0	0	0	0	0	0
# 0	0	0	0	0	0	0	0	0	0	0	0
# 0	0	0	0	0	0	0	0	0	0	0	0
# 0	0	0	0	0	0	0	0	0	0	0	0
# 0	0	0	0	0	0	0	0	0	0	0	0
# 0	0	0	0	0	0	0	0	0	0	0	0

nindivs <- 4
indivs <- paste('ind', seq_len(nindivs)-1, sep="")
bases <- c("A", "C", "G", "T")

angsd_header <- as.vector(t(outer(indivs, bases, paste, sep="_")))

nsnps <- 222

counts <- matrix(rpois(nsnps * length(angsd_header), lambda=3), nrow=nsnps)
colnames(counts) <- angsd_header

As <- 1+4*(0:(nindivs-1))
coverage <- counts[,As]
for (k in 1:3) { coverage <- coverage + counts[,As+k] }
colnames(coverage) <- indivs

countsfile <- "test.counts"
write.table(counts, file=countsfile, row.names=FALSE)

pos <- data.frame( chr=paste0("scaffold_", c(0,0,0, sort(sample.int(6, nsnps-5, replace=TRUE)), 10, 10) ), stringsAsFactors=FALSE)
pos$chr <- factor(pos$chr, levels=unique(pos$chr))
pos$pos <- unlist(lapply(unique(pos$chr), function (x) { cumsum(rpois(sum(pos$chr==x), lambda=5)) }))
pos$totDepth <- rowSums(counts)

posfile <- "test.pos.gz"
poscon <- gzfile(posfile, open="w")
write.table(pos, file=poscon, row.names=FALSE)
close(poscon)

#############
context("Testing counts-to-bin.R")

bincountsfile <- "test.counts.bin"
unlink(bincountsfile)
system(sprintf("../counts-to-bin.R %s %s 77", countsfile, bincountsfile))

expect_true(file.exists(bincountsfile))

# to read lines from the result, NA'ing out counts above the maximum:
get_bincount <- function () {
    bincount <- file(bincountsfile, open="rb")
    attr(bincount,"nindivs") <- nindivs
    attr(bincount,"nbytes") <- 1L
    bincount
}

read_bincounts <- function (bincount, lines) {
    line.coords <- 4*(lines-1) * attr(bincount, "nbytes") * attr(bincount, "nindivs")
    output <- matrix(integer(4 * length(lines) * attr(bincount, "nindivs")), nrow=length(lines))
    for (k in seq_along(lines)) {
        seek(bincount, line.coords[k])
        output[k,] <- readBin( bincount, what=integer(), 
                              n=4*attr(bincount,"nindivs"),
                              size=attr(bincount,"nbytes"), 
                              signed=(attr(bincount,"nbytes")>2) )
    }
    output[output>=256^attr(bincount,"nbytes")-1] <- NA
    return(output)
}

expect_equivalent(as.vector(counts[1,]), as.vector(read_bincounts(get_bincount(), 1)))

check_lines <- c(1,4,5,10,100,101,203)
expect_equivalent(counts[check_lines,], read_bincounts(get_bincount(), check_lines))


##########
context("Testing get-coverage-histogram.R")

max_counts <- 30
truehist <- data.frame(apply( coverage+1L, 2, tabulate, nbins=1+max_counts))

binhist <- read.table(pipe(sprintf("../get-coverage-histogram.R %s 30", bincountsfile)), header=TRUE)

expect_equal(truehist, binhist)

##########
context("Testing get-coverage-by-scaffold.R")

byscaffile <- "test.counts.coverage_by_scaffold"

test_coverage_by_scaffold <- function (min_coverage, max_coverage) {
    unlink(byscaffile)
    blocklines <- 13 #  number of lines to read in at once
    system(sprintf("../get-coverage-by-scaffold.R test.counts.bin %d %d %s %d", min_coverage, max_coverage, byscaffile, blocklines))

    byscaf <- read.table(byscaffile, header=TRUE, stringsAsFactors=FALSE)
    byscaf$scaffold_name <- factor(byscaf$scaffold_name, levels=levels(pos$chr))

    use_these <- (rowSums(counts) > min_coverage) & (rowSums(counts) < max_coverage)
    true_byscaf <- data.frame(scaffold_name=unique(pos$chr))
    true_byscaf$start <- with(subset(pos,use_these), as.vector(tapply(pos, chr, min)))
    true_byscaf$end <- with(subset(pos,use_these), as.vector(tapply(pos, chr, max)))
    true_byscaf$num_sites <- with(subset(pos,use_these), as.vector(table(chr)))
    true_byscaf <- cbind(true_byscaf, 
                         data.frame(do.call(rbind, 
                                        tapply((1:nrow(counts))[use_these], pos$chr[use_these], 
                                               function (k) { colSums(coverage[k,,drop=FALSE]) }))))
    rownames(true_byscaf) <- NULL

    expect_equal(byscaf, true_byscaf)
}

test_coverage_by_scaffold(0,100)
test_coverage_by_scaffold(40,50)
