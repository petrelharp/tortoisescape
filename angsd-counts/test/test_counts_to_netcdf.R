## test counts-to-netcdf.sh

library(ncdf4)

countsfile <- "test.counts.gz"
posfile <- "test.pos.gz"
counts <- as.matrix(read.table(countsfile,header=TRUE))
pos <- read.table(posfile,header=TRUE)

out.truth <- list(
        bases=c("A","C","G","T"),
        indivs=unique(gsub("_[ACGT]","",colnames(counts))),
        chr=pos$chr,
        pos=pos$pos )
out.truth$counts=array(as.numeric(t(counts)),dim=c(4,length(out.truth$indivs),nrow(pos)))

ncfile <- "test.nc"

nin <- nc_open(ncfile)

out.ncdf <- list(
        bases = ncvar_get(nin,"base"),
        indivs = ncvar_get(nin,"indiv"),
        chr = ncvar_get(nin,"chr"),
        pos = ncvar_get(nin,"position"),
        counts = ncvar_get(nin,"count")
    )

nc_close(nin)

for (xn in names(out.truth)) {
    if (!all( out.truth[[xn]] == out.ncdf[[xn]] ) ) { stop(paste(xn, "doesn't match.")) }
}

cat("All good!\n")
