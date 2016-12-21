#!/usr/bin/Rscript

usage <- "
Usage:
   plot-genotypes.R (bincount file) (posfile) (scaffold name) [start:end]
produces a matrix plot of genotypes, with blue for alleles common in the north,
purple for alleles common in the south, red for heterozygotes, etc.

Example:
  plot-genotypes.R 272torts_snp1e6_minmapq20minq30_map2kenro.counts.bin 272torts_snp1e6_minmapq20minq30_map2kenro.pos.gz scaffold_272 
"

arglist <- if (!interactive()) { commandArgs(TRUE) } else { scan(what='char') }
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
    # returns counts where *columns* are sites
    line.coords <- 4*(lines-1)*attr(bincount,"nbytes")*attr(bincount,"nindivs")
    output <- matrix( integer(4*length(lines)*attr(bincount,"nindivs")), ncol=length(lines) )
    for (k in seq_along(lines)) {
        seek(bincount,line.coords[k])
        output[,k] <- readBin( bincount, what=integer(),
                              n=4*attr(bincount,"nindivs"),
                              size=attr(bincount,"nbytes"),
                              signed=(attr(bincount,"nbytes")>2) )
    }
    output[output>=256^attr(bincount,"nbytes")-1] <- NA
    return(output)
}

# read in the matrix of raw counts:
# this is (number of sites) x (4 * number of samples)
counts <- read_bincounts(bincount,pos$line)

# coverages per individual
coverage <- ( counts[1+4*(0:(nindivs-1)),] + counts[2+4*(0:(nindivs-1)),] 
                  + counts[3+4*(0:(nindivs-1)),] + counts[4+4*(0:(nindivs-1)),] )
# coverages per allele across all individuals
totals <- sapply( 1:4, function (k) {
        colSums( counts[k+4*(0:(nindivs-1)),] )
    } )
names(totals) <- c("A","C","G","T")
max.counts <- pmax( totals[,1], totals[,2], totals[,3], totals[,4] )
major.col <- ifelse( totals[,1]==max.counts, 1,
                    ifelse( totals[,2]==max.counts, 2,
                        ifelse( totals[,3]==max.counts, 3, 4 ) ) )
# major allele frequencies by individual
major.coverage <- counts[ cbind( as.vector(outer(4*(0:(nindivs-1)),major.col,"+")), rep(1:ncol(counts),each=nindivs) ) ]
major.freqs <- major.coverage / as.vector(coverage) 
major.freqs[!is.finite(major.freqs)] <- 0
dim(major.freqs) <- dim(major.coverage) <- c(nindivs,ncol(counts))

# color SNPs based on polarization
northern.indivs <- which(pcs$PC1 > 0.04)
southern.indivs <- which(pcs$PC1 < -0.04)

ns.freqs <- rbind( north=colMeans( major.freqs[northern.indivs,] ),
                   south=colMeans( major.freqs[southern.indivs,] ) 
               )

# partition types according to:
#  - if overall frequency is less than 20%, color major allele grey and minor black
#  - if diff in frequency between north and south is more than 0.2,
#    color 'north' type blue and 'south' type purple
#  - otherwise, color allele more common in the north cyan and those in the south plum.

# in this table, first column is color for major allele; second is color for others
snp.type <- ifelse( 
            colMeans(major.freqs)<0.1, 'small',
            ifelse(abs(ns.freqs[1,]-ns.freqs[2,])>=0.2, 'ns', 'other') )

snp.polarization <- sign( ns.freqs['north',]-ns.freqs['south',] )

# matrix of colors corresponding to per-individual major allele
major.colmat <- matrix("",nrow=nrow(major.freqs),ncol=ncol(major.freqs))
major.colmat[,snp.type=='small'] <- ifelse(major.freqs[,snp.type=='small']>0.5,'grey','black')
major.colmat[,snp.type=='ns'] <- ifelse(((major.freqs[,snp.type=='ns']-0.5)*snp.polarization[snp.type=='ns'])>=0,'blue','purple')
major.colmat[,snp.type=='other'] <- ifelse(((major.freqs[,snp.type=='other']-0.5)*snp.polarization[snp.type=='other'])>=0,'cyan','plum')
major.colmat[coverage==0] <- NA

# matrix of colors corresponding to per-individual minor allele if present
minor.colmat <- matrix("",nrow=nrow(major.freqs),ncol=ncol(major.freqs))
minor.colmat[,snp.type=='small'] <- ifelse((1-major.freqs)[,snp.type=='small']>0.5,'grey','black')
minor.colmat[,snp.type=='ns'] <- ifelse((0.5-major.freqs[,snp.type=='ns'])*(ns.freqs[1,snp.type=='ns']-ns.freqs[2,snp.type=='ns'])>=0,'blue','purple')
minor.colmat[,snp.type=='other'] <- ifelse((0.5-major.freqs[,snp.type=='other'])*(ns.freqs[1,snp.type=='other']-ns.freqs[2,snp.type=='other'])>=0,'cyan','plum')
minor.colmat[(coverage-major.coverage)==0] <- NA

# plotting order
plot.ord <- order(pcs$PC1)

# png(file="mt-haplotypes.png",width=5*288,height=2.5*288,res=288,pointsize=8)
opar <- par(mar=c(1,1,8,1)+.1)
plot( row(major.colmat[plot.ord,]), col(major.colmat[plot.ord,]), 
     main=sprintf("haplotypes on %s",scaffold),
        col=adjustcolor(major.colmat[plot.ord,],0.75),
        bg=adjustcolor(minor.colmat[plot.ord,],0.75),
        pch=21, cex=0.5,
        xlab='', xaxt='n', yaxt='n', ylab='' )
axis(3, at=seq_along(tort.ids), labels=gsub(" (sheared2)","",tort.ids[plot.ord],fixed=TRUE), cex=0.25, las=3, cex.axis=0.25 )
par(opar)
