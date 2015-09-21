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
geog <- read.csv("../tortoisescape/tort_272_info/geog_distance.csv")
etorts <- unique(c(levels(geog$etort1),levels(geog$etort2)))
netorts <- as.numeric(gsub("[^0-9]","",gsub(" .*","",etorts)))
etorts <- etorts[order(netorts)]
ntorts <- length(etorts)

ibd$etort1 <- factor( etorts[match(ibd$id1,ids)], levels=etorts )
ibd$etort2 <- factor( etorts[match(ibd$id2,ids)], levels=etorts )

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

pdf(file="ibd-by-distance.pdf",width=10,height=8,pointsize=10)
matplot(x=dist.mids, y=len.dist, type='l', log='y', ylab="mean number of blocks", xlab="geographical distance", col=rainbow(length(len.mids)), lty=1 )
legend("topright",lty=1,col=rainbow(length(len.mids)), legend=paste(levels(ibd$length.bin),"bp"))
dev.off()
