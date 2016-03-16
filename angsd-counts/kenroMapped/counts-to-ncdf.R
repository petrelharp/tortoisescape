#!env Rscript

# convert the .counts file to a netcdf file to allow quick access of individual sites

library(ncdf4)

tortdir <- gsub("tortoisescape.*","tortoisescape",getwd())

datadir <- file.path(tortdir,"angsd-counts/kenroMapped")
countfile <- file.path(datadir,"272torts_snp1e6_minmapq20minq30_map2kenro.counts.gz")

# output file:
ncfile <- gsub(".gz$",".nc",countfile)
cat("Writing to ", ncfile, "\n")

count.con <- pipe(paste("zcat",countfile),open="r")
count.header <- scan(count.con,nlines=1,what="char")

pos <- read.table(file.path(datadir,"272torts_snp1e6_minmapq20minq30_map2kenro.pos.gz"),header=TRUE)

nindivs <- length(count.header)/4

# note this will be TRANSPOSED to whats in the counts file
ypos <- seq_along(count.header)/4
dim1 <- ncdim_def('individual', 'samples', as.double(ypos))

# make fake coordinates out of position
maxpos <- tapply( pos$pos, pos$chr, max, na.rm=TRUE )
cumpos <- c(0,cumsum(as.numeric(maxpos)))
xpos <- cumpos[as.numeric(pos$chr)] + pos$pos/maxpos[as.numeric(pos$chr)]
dim2 <- ncdim_def('position', 'snps', as.double(xpos))

varz <- ncvar_def( name='Allele counts',
                   units='number', 
                   dim=list(dim1, dim2), 
                   missval=NULL,
                   longname = paste("allele counts from",countfile) )
 
outnc <- nc_create(ncfile, varz, force_v4=TRUE)
 
blocksize <- 1e4

ncatt_put(outnc, 0, "colnames", count.header)
ncatt_put(outnc, 0, "chr", pos$chr)

nlines <- 0
while (length(
              counts <- scan(count.con,nlines=blocksize)
              )>0) {
    dim(counts) <- c(4*nindivs,length(counts)/(4*nindivs))
    # rownames(counts) <- count.header
    ncvar_put(outnc, varz, counts, start=c(1,1+nlines), count=c(-1,blocksize))
    nlines <- nlines+ncol(counts)
}
nc_close(outnc)


if (FALSE) {
    ## TESTING: ncdf writes in column order, same as R
    xx <- matrix( 1:120, nrow=10 )
    library(ncdf4)
    dim1 <- ncdim_def("x","m",1:nrow(xx))
    dim2 <- ncdim_def("y","n",1:ncol(xx))
    varz <- ncvar_def(name="test",units='unit',dim=list(dim1,dim2),missval=NULL)
    outnc <- nc_create("test.nc",varz,force_v4=TRUE)
    ncvar_put(outnc, varz, xx[1:60],
              start=c(1,1), count=c(-1,6))
    ncvar_put(outnc, varz, xx[60+(1:60)],
              start=c(1,7), count=c(-1,6))
    nc_close(outnc)
    ##
    ncin <- nc_open("test.nc")
    yy <- ncvar_get(ncin,"test")
    stopifnot(all.equal(xx,yy))
}

