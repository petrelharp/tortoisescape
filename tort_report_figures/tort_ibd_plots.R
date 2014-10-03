################################################################
################################################################
#	make various IBD-themed plots for tortoises
################################################################
################################################################

################################
#	some preliminaries
################################

if(file.exists("~/Desktop/Dropbox/tortoisescape")){
	setwd("~/Desktop/Dropbox/tortoisescape")
}

# get the metadata
torts <- read.csv("1st_180_torts.csv",header=TRUE,stringsAsFactors=FALSE)
nind <- nrow(torts)

# and the pairwise geographic distance matrix
tort.dist.table <- read.table("1st180_pairwise_distances_sorted_redundancy_removed.txt",header=TRUE)
tort.dists <- numeric(nind^2); dim(tort.dists) <- c(nind,nind)
tort.dists[ cbind( match(tort.dist.table$etort1,torts$EM_Tort_ID), 
					match(tort.dist.table$etort2,torts$EM_Tort_ID) ) ] <- tort.dist.table$DISTANCE
tort.dists <- tort.dists + t(tort.dists)


# then grab the pairwise pi matrix
tort.pwp.vals <- scan("pairwisePi/alleleCounts_1millionloci.pwp")
tort.pwp <- numeric(nind^2)
dim(tort.pwp) <- c(nind,nind)
	tort.pwp[upper.tri(tort.pwp,diag=TRUE)] <- tort.pwp.vals
	tort.pwp[lower.tri(tort.pwp,diag=FALSE)] <- t(tort.pwp)[lower.tri(tort.pwp,diag=FALSE)]
dimnames(tort.pwp) <- list( torts$EM_Tort_ID, torts$EM_Tort_ID )

# do PCA on tortoise genetic covariance
covmat <- as.matrix(read.table("covmat/alleleCounts500kLoci-covmat.txt"))
pmat <- diag(nind) - 1/nind
eig.covmat <- eigen(pmat %*% covmat %*% pmat)

# use the PCA to make plotting color & symbol options
tort.plotting.colors.discrete <- rep("red",nind)
	tort.plotting.colors.discrete[eig.covmat$vectors[,1] > 0] <- "purple"
tort.plotting.colors.continuous <- rev(rainbow(nind,start=4/6,end=6/6))[as.numeric(cut(eig.covmat$vectors[,1],nind))]
tort.plot.mat.col.2bins <- matrix("gray",nrow=nind,ncol=nind)
tort.plot.mat.col.3bins <- matrix("gray",nrow=nind,ncol=nind)
	for(i in 1:nind){
		for(j in 1:nind){
			if(tort.plotting.colors.discrete[i] == tort.plotting.colors.discrete[j]){
				if(tort.plotting.colors.discrete[i] == "red"){
					tort.plot.mat.col.3bins[i,j] <- "red"
				} else if(tort.plotting.colors.discrete[i] == "purple"){
					tort.plot.mat.col.3bins[i,j] <- "purple"
				}
			} else if(tort.plotting.colors.discrete[i] != tort.plotting.colors.discrete[j]){
				tort.plot.mat.col.2bins[i,j] <- "green"
				tort.plot.mat.col.3bins[i,j] <- "green"
			}
		}
	}

################################
#	FIGURES
################################

index.mat <- upper.tri(tort.dists,diag=FALSE)

################
#	IBD colored by
#		N-N,N-S,S-S
################

png(file="tort_report_figures/IBD_NN_SS_NS_plot.png",res=200,width=6*200,height=5*200)
	plot(tort.dists[index.mat],tort.pwp[index.mat],
			col=tort.plot.mat.col.3bins[index.mat],
			pch=20,cex=0.3,xlab="pairwise geographic distance",ylab="pairwise sequence divergence",main = "Isolation by Distance")
	for(i in 1:length(unique(c(tort.plot.mat.col.3bins)))){
		use.these <- which(tort.plot.mat.col.3bins == unique(c(tort.plot.mat.col.3bins))[i],arr.ind=TRUE)
		lines(lowess(tort.pwp[use.these] ~ tort.dists[use.these]),col=unique(c(tort.plot.mat.col.3bins))[i],lwd=2)
	}
	legend(x="bottomright",lty=1,col=unique(c(tort.plot.mat.col.3bins)),legend=c("North - North","North - South","South - South"))
	#points(diag(tort.dists),diag(tort.pwp),pch=20,cex=0.5)
dev.off()

################
#	IBD colored by
#		same-diff
################

png(file="tort_report_figures/IBD_same_diff_plot.png",res=200,width=6*200,height=5*200)
	plot(tort.dists[index.mat],tort.pwp[index.mat],
			col=tort.plot.mat.col.2bins[index.mat],
			pch=20,cex=0.3,xlab="pairwise geographic distance",ylab="pairwise sequence divergence",main = "Isolation by Distance")
	for(i in 1:length(unique(c(tort.plot.mat.col.2bins)))){
		use.these <- which(tort.plot.mat.col.2bins == unique(c(tort.plot.mat.col.2bins))[i],arr.ind=TRUE)
		lines(lowess(tort.pwp[use.these] ~ tort.dists[use.these]),col=unique(c(tort.plot.mat.col.2bins))[i],lwd=2)
	}
	legend(x="bottomright",lty=1,col=unique(c(tort.plot.mat.col.2bins)),legend=c("North - North and South-South","North - South"))
dev.off()


################
#	IBD colored by
#		distance on PC1
################

pc1.dist.mat <- fields::rdist(eig.covmat$vectors[,1])
pc1.dist.mat.cols <- heat.colors(length(pc1.dist.mat))[as.numeric(cut(pc1.dist.mat,length(pc1.dist.mat)))]
png(file="tort_report_figures/IBD_PC1dist_plot.png",res=200,width=6*200,height=5*200)
	#quartz(width=6,height=5)
	par(bg="gray")
	plot(tort.dists[index.mat],tort.pwp[index.mat],
			col=pc1.dist.mat.cols[index.mat],
			pch=20,cex=0.3,xlab="pairwise geographic distance",ylab="pairwise sequence divergence",main = "Isolation by Distance")
				leg.x <- c(1:nind)
				leg.y <- sort(pc1.dist.mat[index.mat])
				leg.z <- outer(leg.x,leg.y)
				leg.lab.vals <- c(round(min(pc1.dist.mat[index.mat]),2),round(max(pc1.dist.mat[index.mat]),2))
				fields::image.plot(leg.x,leg.y,leg.z,add=TRUE,
						nlevel=nind,
						horizontal=FALSE,graphics.reset=TRUE,
						legend.only=TRUE,col=heat.colors(length(pc1.dist.mat)),
						smallplot=c(0.815,0.835,0.25,0.5),
						axis.args=list(at=c(0.03,35.4),labels=leg.lab.vals))
				mtext(text="Distance on\n genetic PC1",side=1,padj=-2.7,adj=0.78)
dev.off()


################
#	correlation in pairwise pi
#		against distance
################

pi.cor <- cor(tort.pwp)

png(file="tort_report_figures/corr_IBD_NN_SS_NS_plot.png",res=200,width=6*200,height=5*200)
	#quartz(width=6,height=5)
#	par(bg="gray")
	plot(tort.dists[index.mat],pi.cor[index.mat],
			col=tort.plot.mat.col.3bins[index.mat],
			pch=20,cex=0.3,xlab="pairwise geographic distance",
			ylab="pairwise sequence divergence",main = "Correlation in Isolation by Distance")
	for(i in 1:length(unique(c(tort.plot.mat.col.3bins)))){
		use.these <- which(tort.plot.mat.col.3bins == unique(c(tort.plot.mat.col.3bins))[i],arr.ind=TRUE)
		lines(lowess(pi.cor[use.these] ~ tort.dists[use.these]),col=unique(c(tort.plot.mat.col.3bins))[i],lwd=2)
	}
	legend(x="topright",lty=1,col=unique(c(tort.plot.mat.col.3bins)),
			legend=c("North - North","North - South","South - South"),cex=0.7)
dev.off()

png(file="tort_report_figures/corr_IBD_same_diff_plot.png",res=200,width=6*200,height=5*200)
	#quartz(width=6,height=5)
#	par(bg="gray")
	plot(tort.dists[index.mat],pi.cor[index.mat],
			col=tort.plot.mat.col.2bins[index.mat],
			pch=20,cex=0.3,xlab="pairwise geographic distance",
			ylab="pairwise sequence divergence",main = "Correlation in Isolation by Distance")
	for(i in 1:length(unique(c(tort.plot.mat.col.2bins)))){
		use.these <- which(tort.plot.mat.col.2bins == unique(c(tort.plot.mat.col.2bins))[i],arr.ind=TRUE)
		lines(lowess(pi.cor[use.these] ~ tort.dists[use.these]),col=unique(c(tort.plot.mat.col.2bins))[i],lwd=2)
	}
	legend(x="topright",lty=1,col=unique(c(tort.plot.mat.col.2bins)),legend=c("North - North and South-South","North - South"),cex=0.7)
dev.off()

png(file="tort_report_figures/corr_IBD_PC1distcol_plot.png",res=200,width=6*200,height=5*200)
	#quartz(width=6,height=5)
	par(bg="gray")
	plot(tort.dists[index.mat],pi.cor[index.mat],
			col=pc1.dist.mat.cols[index.mat],
			pch=20,cex=0.3,xlab="pairwise geographic distance",
			ylab="pairwise sequence divergence",main = "Correlation in Isolation by Distance")
				leg.x <- c(1:nind)
				leg.y <- sort(pc1.dist.mat[index.mat])
				leg.z <- outer(leg.x,leg.y)
				leg.lab.vals <- c(round(min(pc1.dist.mat[index.mat]),2),round(max(pc1.dist.mat[index.mat]),2))
				fields::image.plot(leg.x,leg.y,leg.z,add=TRUE,
						nlevel=nind,
						horizontal=FALSE,graphics.reset=TRUE,
						legend.only=TRUE,col=heat.colors(length(pc1.dist.mat)),
						smallplot=c(0.815,0.835,0.5,0.75),
						axis.args=list(at=c(0.03,35.4),labels=leg.lab.vals))
				mtext(text="Distance on\n genetic PC1",side=1,padj=-6.7,adj=0.78)
dev.off()