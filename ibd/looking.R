#   1)  First sample identifier
#   2)  First sample haplotype index (1 or 2)
#   3)  Second sample identifier
#   4)  Second sample haplotype index (1 or 2)
#   5)  Chromosome
#   6)  Starting genomic position (inclusive)
#   7)  Ending genomic position (inclusive)
#   8)  LOD score (larger values indicate greater evidence for IBD)

ibd <- read.table("beagleAllSNPs.ibd",stringsAsFactors=FALSE)
names(ibd) <- c("id1","hap1","id2","hap2","chrom","start","end","LOD")
ids <- unique(c(ibd$id1,ibd$id2))
nids <- as.numeric( gsub("ind","",ids) )
ids <- ids[order(nids)]
geog <- read.csv("../tort_272_info/geog_distance.csv")
pcs <- read.csv("../tort_272_info/pcs.csv")
etorts <- unique(c(levels(geog$etort1),levels(geog$etort2)))
netorts <- as.numeric(gsub("[^0-9]","",gsub(" .*","",etorts)))
etorts <- etorts[order(netorts)]
ntorts <- length(etorts)
chroms <- unique(ibd$chrom)

ibd$etort1 <- factor( etorts[match(ibd$id1,ids)], levels=etorts )
ibd$etort2 <- factor( etorts[match(ibd$id2,ids)], levels=etorts )
ibd$chrom <- factor( ibd$chrom, levels=chroms )

chrlens <- with(ibd,tapply(end,chrom,max)-tapply(start,chrom,min))

geog.mat <- matrix( nrow=ntorts, ncol=ntorts )
rownames(geog.mat) <- colnames(geog.mat) <- etorts
geog.mat[ cbind( match(geog$etort1,etorts), match(geog$etort2,etorts) ) ] <- geog$distance
geog.mat[is.na(geog.mat)] <- t(geog.mat)[is.na(geog.mat)]
ut <- upper.tri(geog.mat,diag=FALSE)

count.ibds <- function (subs=TRUE) {
    tab <- with( subset(ibd,subs), table(etort1,etort2) )
    tab <- (tab+t(tab)-diag(diag(tab)))
    tab[ut]
}
all.tab <- count.ibds()
stopifnot(all(rownames(all.tab)==etorts))


long.tab <- count.ibds(ibd$LOD>5 & (ibd$end-ibd$start>1e4))

nibd <- data.frame(
        etort1=etorts[row(geog.mat)[ut]],
        etort2=etorts[col(geog.mat)[ut]],
        dist=geog.mat[ut]
        )

len.breaks <- c(0,exp(seq(log(1000),log(2e4),length.out=10)),1e5)
len.mids <- len.breaks[-1] - diff(len.breaks)/2
len.mids[length(len.mids)] <- rev(len.breaks)[2]
ibd$length.bin <- cut( ibd$end-ibd$start,breaks=len.breaks )
length.bin.names <- sapply(1:(length(len.breaks)-1), function (k) {
            sprintf("%0.1f-%0.1f kB", len.breaks[k]/1000, len.breaks[k+1]/1000)
        } )

for (k in 1:nlevels(ibd$length.bin)) {
    nibd[[levels(ibd$length.bin)[k]]] <- count.ibds( as.numeric(ibd$length.bin)==k )
}

dist.breaks <- seq(0,max(geog.mat),length.out=20)
dist.mids <- dist.breaks[-1]-diff(dist.breaks)/2
nibd$dist.bin <- cut( nibd$dist, breaks=dist.breaks )

len.dist <- sapply( 1:nlevels(ibd$length.bin), function (k) {
                tapply( nibd[[levels(ibd$length.bin)[k]]], nibd$dist.bin, mean )
        } )
colnames(len.dist) <- levels(ibd$length.bin)

pdf(file="ibd-by-distance.pdf",width=10,height=8,pointsize=10)
matplot( x=dist.mids, y=len.dist, type='l', log='y', ylab="mean number of blocks", xlab="geographical distance", col=rainbow(length(len.mids)+5), lty=1 )
legend("topright",lty=1,col=rainbow(length(len.mids)+5), legend=paste(levels(ibd$length.bin),"bp"))
dev.off()

pdf(file="ibd-by-distance-for-talk.pdf",width=5,height=3.5,pointsize=10)
matplot( x=dist.mids/1000, y=len.dist[,len.breaks[-1]>4e3], type='l', log='y', xlim=c(0,600),
        ylab="mean number of blocks", xlab="geographical distance (km)", col=rainbow(length(len.mids)+5), lty=1,
       main="IBD blocks per pair" )
legend("bottomright",lty=1,col=rainbow(length(len.mids)+5), legend=length.bin.names[len.breaks[-1]>4e3], bg="white", cex=0.75)
dev.off()

## length spectrum by distance

nibd <- data.frame(
        etort1=etorts[row(geog.mat)[ut]],
        etort2=etorts[col(geog.mat)[ut]],
        dist=geog.mat[ut]
        )


len.breaks <- c(0,exp(seq(log(1000),log(2e4),length.out=50)),1e5)
len.mids <- len.breaks[-1] - diff(len.breaks)/2
len.mids[length(len.mids)] <- rev(len.breaks)[2]
ibd$length.bin <- cut( ibd$end-ibd$start,breaks=len.breaks )

for (k in 1:nlevels(ibd$length.bin)) {
    nibd[[levels(ibd$length.bin)[k]]] <- count.ibds( as.numeric(ibd$length.bin)==k )
}

dist.breaks <- seq(0,max(geog.mat),length.out=20)
dist.mids <- dist.breaks[-1]-diff(dist.breaks)/2
nibd$dist.bin <- cut( nibd$dist, breaks=dist.breaks )

len.dist <- sapply( 1:nlevels(ibd$length.bin), function (k) {
                tapply( nibd[[levels(ibd$length.bin)[k]]], nibd$dist.bin, mean )
        } )
colnames(len.dist) <- levels(ibd$length.bin)

matplot(len.mids, t(len.dist),type='l',log='y')

### by north/south

ns.labels <- ifelse( pcs$PC1[match(etorts,pcs$etort)]>0, "N", "S" )
nibd$ns <- paste( ns.labels[match(nibd$etort1,etorts)], ns.labels[match(nibd$etort2,etorts)], sep="")
ns.spectra <- sapply( 1:nlevels(ibd$length.bin), function (k) {
                tapply( nibd[[levels(ibd$length.bin)[k]]], nibd$ns, mean )
        } )
colnames(ns.spectra) <- levels(ibd$length.bin)

pdf(file="IBD-spectrum-by-NS.pdf",width=4,height=4,pointsize=10)
matplot(len.mids, t(ns.spectra),type='l',log='y',lty=1:4,col=1:4,xlab="IBD length (bp)", ylab="proportion")
legend("topright",legend=rownames(ns.spectra),lty=1:4,col=1:4)
dev.off()

####

require(Matrix)

indpair <- function (etort1,etort2) {
    # return a unique index for an unordered pair of etorts
    ii <- matrix(1:length(etorts)^2,nrow=length(etorts))
    ut <- upper.tri(ii,diag=TRUE)
    ii[ut] <- (1:sum(ut))
    ii[lower.tri(ii,diag=FALSE)] <- t(ii)[lower.tri(ii,diag=FALSE)]
    return( ii[ cbind(match(etort1,etorts),match(etort2,etorts)) ] )
}

etort.pairs <- outer(etorts,etorts,paste,sep="-")
etort.pairs <- etort.pairs[upper.tri(etort.pairs,diag=TRUE)]
cii <- with(ibd, sparseMatrix( i=indpair(etort1,etort2), j=as.numeric(chrom), x=len, dims=c(length(etort.pairs),length(chroms)) ) )

contig.means <- colMeans(cii>0)
pair.means <- rowMeans(cii>0)

# (X-A)^T (X-A) = X^T X - X^T A - A^T X + A^T A
#  X^T (u 1^T) = X^T u 1^T
#  X^T (1 v^T) = X^T 1 v^T

contig.cov <- crossprod(cii>0) - outer(colMeans(cii>0),colMeans(cii>0))
contig.cor <- contig.cov / (nrow(cii)*sqrt( outer(colMeans(cii>0),colMeans(cii>0)) ))
diag(contig.cor) <- 0

perfects <- which(contig.cor>0.9 & row(contig.cor)>col(contig.cor))
cbind( row(contig.cor)[perfects], col(contig.cor)[perfects] )

