################################################################
################################################################
#	make IBD&E plots for tortoises
################################################################
################################################################

################################
#	some preliminaries
################################

if(file.exists("~/Desktop/Dropbox/tortoisescape/tort_report_figures")){
	setwd("~/Desktop/Dropbox/tortoisescape/tort_report_figures")
}

# get the metadata
torts <- read.csv("../1st_180_torts.csv",header=TRUE,stringsAsFactors=FALSE)
nind <- nrow(torts)

# and the pairwise geographic distance matrix
tort.dist.table <- read.table("../1st180_pairwise_distances_sorted_redundancy_removed.txt",header=TRUE)
tort.dists <- numeric(nind^2); dim(tort.dists) <- c(nind,nind)
tort.dists[ cbind( match(tort.dist.table$etort1,torts$EM_Tort_ID), 
					match(tort.dist.table$etort2,torts$EM_Tort_ID) ) ] <- tort.dist.table$DISTANCE
tort.dists <- tort.dists + t(tort.dists)


# then grab the pairwise pi matrix
tort.pwp.vals <- scan("../pairwisePi/alleleCounts_1millionloci.pwp")
tort.pwp <- numeric(nind^2)
dim(tort.pwp) <- c(nind,nind)
	tort.pwp[upper.tri(tort.pwp,diag=TRUE)] <- tort.pwp.vals
	tort.pwp[lower.tri(tort.pwp,diag=FALSE)] <- t(tort.pwp)[lower.tri(tort.pwp,diag=FALSE)]
dimnames(tort.pwp) <- list( torts$EM_Tort_ID, torts$EM_Tort_ID )

# now, go through each edist .tsv file
#	and make an IBD plot colored by edist

################################
# some quick functions to help out
################################

# makes a symmetric matrix out from a
#	list of pairwise e-distances

matrixify.edist.tsv <- function(tsv.file,nind,EM_Tort_ID){
	edist.vals <- read.table(tsv.file,header=TRUE)
	edist.mat <- numeric(nind^2)
	dim(edist.mat) <- c(nind,nind)
	edist.mat[cbind( match(edist.vals$tort1,EM_Tort_ID),
					match(edist.vals$tort2,EM_Tort_ID) ) ] <- edist.vals$edist
	return(edist.mat)
}

# generates a vector of plotting colors 
#	from the e-distances

make.edist.plotting.colors <- function(edist.mat,index.mat){
	n.cols <- length(edist.mat[index.mat])
	edist.plotting.cols <- adjustcolor(diverge_hcl(n.cols,c=100),.5)[as.numeric(cut(order(edist.mat[index.mat]),n.cols))]   # [as.numeric(cut(edist.mat[index.mat],n.cols))]
	return(edist.plotting.cols)
}


################################
# go through each .tsv and make
#	and IBD&E .png
################################

if (file.exists("~/Desktop/10x")){
	tsvdir <- "~/Desktop/10x"
} else {
    tsvdir <- "../geodist/10x"
}

outdir <- "ibde.pngs/"

index.mat <- upper.tri(tort.pwp,diag=FALSE)
tsv.files <- list.files(tsvdir,pattern="*.tsv",full.names=TRUE)
tsv.file.names <- gsub(".tsv","",tsv.files)

tort.pwp.ut <- tort.pwp[index.mat]
tort.dists.ut <- tort.dists[index.mat]

if(!file.exists(outdir)){
	dir.create(outdir)
}
for(i in 1:length(tsv.files)){
	cat(i,"\t")
	edist.mat <- matrixify.edist.tsv(tsv.files[i],nind,torts$EM_Tort_ID)
	edist.plotting.cols <- make.edist.plotting.colors(edist.mat,index.mat)
	png(file=paste(outdir,basename(tsv.file.names[i]),"ibde.png",sep=""),res=200,width=8*200,height=5*200)
		plot(tort.dists.ut,tort.pwp.ut,col=edist.plotting.cols,pch=20,cex=0.5,
				xlab="pairwise geo distance",
				ylab="pairwise pi",
				main=paste("IBD&E: ",basename(tsv.file.names[i]),sep=""))
        edist.qs <- quantile(edist.mat[index.mat],c(0,.25,.5,.75,1),na.rm=TRUE)
        qcols <- rev(diverge_hcl(4,c=100))
        for (k in 1:4) {
            usethese <- (!is.na(edist.mat[index.mat])) & ( edist.mat[index.mat] >= edist.qs[k] ) & ( edist.mat[index.mat] < edist.qs[k+1] )
            if (any(usethese)) { lines( lowess( tort.dists.ut[usethese], tort.pwp.ut[usethese], f=.2 ), col=qcols[k], lwd=2 ) }
        }
        legend("bottomright", lwd=2, col=qcols, legend=paste(c("first","second","third","fourth"),"quartile"), title=basename(tsv.file.names[i]) )
	dev.off()
}



if(FALSE){
	require(rgl)
	index.mat <- upper.tri(tort.dists,diag=FALSE)
	plot3d(tort.dists[index.mat],slope.edist[index.mat],tort.pwp[index.mat],col=plotting.colors)
	play3d(spin3d())


slope.edist.vals <- read.table("~/Desktop/slope_30-mean-edist.tsv",header=TRUE)
slope.edist <- numeric(nind^2)
dim(slope.edist) <- c(nind,nind)
slope.edist[cbind( match(slope.edist.vals$tort1,torts$EM_Tort_ID),
					match(slope.edist.vals$tort2,torts$EM_Tort_ID) ) ] <- slope.edist.vals[,3]

plotting.colors <- heat.colors(
						length(slope.edist[index.mat]))[as.numeric(
								cut(slope.edist[index.mat],
										length(slope.edist[index.mat])))]

plot(tort.dists[index.mat], tort.pwp[index.mat],col=1,pch=19,cex=0.25)
	points(tort.dists[index.mat], tort.pwp[index.mat],col=plotting.colors,pch=19,cex=0.25)

#tortoise 58 is out of bounds!!!
tort.pwp.e <- tort.pwp[-match("etort-58",row.names(tort.pwp)),
						-match("etort-58",row.names(tort.pwp))]
tort.dists.e <- tort.dists[-match("etort-58",row.names(tort.pwp)),
							-match("etort-58",row.names(tort.pwp))]
}





