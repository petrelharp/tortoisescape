#!env Rscript

library(ncdf4)

datadir <- file.path(tortdir,"angsd-counts/kenroMapped")
countfile <- file.path(datadir,"272torts_snp1e6_minmapq20minq30_map2kenro.counts.gz")

count.con <- pipe(paste("zcat",countfile),open="r")
count.header <- scan(count.con,nlines=1,what="char")
pos <- read.table(file.path(datadir,"272torts_snp1e6_minmapq20minq30_map2kenro.pos.gz"),header=TRUE)

# make fake coordinates out of position
maxpos <- tapply( pos$pos, pos$chr, max, na.rm=TRUE )
cumpos <- cumsum(maxpos)
xpos <- cumpos[as.numeric(pos$chr)] + pos$pos/maxpos[as.numeric(pos$chr)]
ypos <- seq_along(count.header)/4

nindivs <- length(count.header)/4

dim1 <- ncdim_def('position', 'snps', as.double(xpos))
dim2 <- ncdim_def('individual', 'samples', as.double(ypos))

varz <- ncvar_def( name='Allele counts',
                   units='number', 
                   dim=list(dim1, dim2), 
                   missval=NULL,
                   longname = paste("allele counts from",countfile) )
 
outnc <- nc_create(gsub(".gz$",".nc",countfile) varz, force_v4 = TRUE)
 
blocksize <- 1e4

nlines <- 1
while (length(counts <- scan(count.con,nlines=blocksize))>0) {
    dim(counts) <- c(4*nindivs,length(counts)/(4*nindivs))
    # rownames(counts) <- count.header
    ncvar_put(outnc, varz, counts, start=c(1,nlines), count=c(-1,blocksize))
    nlines <- nlines+ncol(counts)
}

nc_close(outnc)
