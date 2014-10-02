################################################################
################################################################
#	visualize tortoise spacemix output
################################################################
################################################################

# specify which analyses to look at
#	and load relevant files
setwd("spacemix/locations")
setwd("locations_2")
load(list.files(pattern="MCMC.output"))
load(list.files(pattern="tort180_spacemix_dataset.Robj"))

# also specify burnin & thinning
burnin <- 0.2
n.mcmc.samples <- 200
k <- last.params$k

################################
#	functions
################################

procrusteez <- function(obs.locs,target.locs,k,source.locs = NULL,option){
	require(vegan)
	proc.loc <- procrustes(obs.locs,target.locs,scale=TRUE)
	if(option==1){
		proc.pop.loc <- proc.loc$scale * target.locs %*% proc.loc$rotation + matrix(proc.loc$translation,nrow=k,ncol=2,byrow=TRUE)
	} else if(option==2){
		proc.pop.loc <- proc.loc$scale * source.locs %*% proc.loc$rotation + matrix(proc.loc$translation,nrow=k,ncol=2,byrow=TRUE)
	}
	return(proc.pop.loc)	
}

fade.admixture.source.points <- function(pop.cols,admix.proportions){
	faded.colors <- numeric(length(pop.cols))
	for(i in 1:length(pop.cols)){
		faded.colors[i] <- adjustcolor(pop.cols[i],admix.proportions[i])
	}
	return(faded.colors)
}

################################
#	get metadata for plotting
################################

# get the tortoise names from Evan's sample numbers
tort.names <- unlist(lapply(strsplit(row.names(tort.coords),"-"),"[",2))

# mean center the covariance matrix and do quick pca on genetic data
pmat <- diag(k) - 1/k
eig.covmat <- eigen(pmat %*% sample.covariance %*% pmat)

# color tortoises by same 
tort.plot.cols <- rep("purple",k)
	tort.plot.cols[eig.covmat $vectors[,1] > 0] <- "blue"

################################
#	procrustes transform estimated 
#		coordinates around observed
# 	and weight fading of colors by 
#		admix proportion
################################
x <- seq((length(which(Prob!=0))*burnin)+1,length(which(Prob!=0)),length.out=n.mcmc.samples)
procrustes.target.locations <- vector("list",length(x))
procrustes.source.locations <- vector("list",length(x))
point.tort.col <- matrix(0,nrow=k,ncol=length(x))

	for(i in 1:length(x)){
		procrustes.target.locations[[i]] <- procrusteez(obs.locs = tort.coords,
											target.locs = population.coordinates[[x[i]]][1:k,],
											k = k,
											option = 1)
		procrustes.source.locations[[i]] <- procrusteez(obs.locs = tort.coords,
											target.locs = population.coordinates[[x[i]]][1:k,],
											k = k,
											source.locs = population.coordinates[[x[i]]][(k+1):(2*k),],
											option = 2)
		point.tort.col[,i] <- fade.admixture.source.points(tort.plot.cols,admix.proportions[,x[i]])
	}

x.min <- min(unlist(lapply(procrustes.target.locations,"[",,1)))
x.max <- max(unlist(lapply(procrustes.target.locations,"[",,1)))
y.min <- min(unlist(lapply(procrustes.target.locations,"[",,2)))
y.max <- max(unlist(lapply(procrustes.target.locations,"[",,2)))

plot(tort.coords,type='n',xlim=c(x.min,x.max),ylim=c(y.min,y.max))
	points(tort.coords,pch=15,col=tort.plot.cols)
	points(procrustes.target.locations[[length(x)]][1:k,],col=tort.plot.cols,pch=20)
	arrows(x0 = procrustes.target.locations[[length(x)]][1:k,1],
			y0 = procrustes.target.locations[[length(x)]][1:k,2], 
			x1 = tort.coords[,1],
			y1 = tort.coords[,2],
			col="green",
			length=0.1,
			lty=1,
			lwd=0.7)


plot(tort.coords,type='n',
		xlim=c(x.min,x.max),
		ylim=c(y.min,y.max),
	for(i in 1:length(x)){
		points(procrustes.target.locations[[i]],pch=20,col=adjustcolor(tort.plot.cols,0.05))
	}


tort.plot.cols2 <- rev(rainbow(k,start=4/6,end=6/6))[as.numeric(cut(tort.coords[,2],k))]

png(file="tort_location_spacemix.png",res=200,height=5*200,width=15*200)
	par(mfrow=c(1,3))
	plot(tort.coords,type='n',
			xlab="Long",ylab="Lat",
			main="Tortoise Locations")
		text(tort.coords,labels=tort.names,col=tort.plot.cols2)
	plot(last.params$population.coordinates[1:k,],type='n',
			ylab="'lat'",xlab="'long'",
			main="Tortoise SpaceMix Locations")
		text(last.params$population.coordinates[1:k,],
				labels= tort.names,
				col= tort.plot.cols2)
	plot(last.params$population.coordinates[1:k,],type='n',
			ylab="'lat'",xlab="'long'",
			main="Tortoise SpaceMix Locations")
		text(last.params$population.coordinates[1:k,],
				labels= tort.names,
				col= tort.plot.cols)
dev.off()


point.tort.ad.cols <- tort.plot.cols2
for(i in 1:k){
	point.tort.ad.cols[i] <- adjustcolor(tort.plot.cols2[i],last.params$admix.proportions[i])
}

load("~/Desktop/Dropbox/tortoisescape/spacemix/admixture_locations/admixture_locations_2/tort_admixture_locations_2space_MCMC_output1.Robj")
load(list.files(pattern="tort180_spacemix_dataset.Robj"))
burnin <- 0.2
n.mcmc.samples <- 200
k <- last.params$k

tort.plot.cols2 <- rev(rainbow(k,start=4/6,end=6/6))[as.numeric(cut(tort.coords[,2],k))]
point.tort.ad.cols <- tort.plot.cols2
for(i in 1:k){
	point.tort.ad.cols[i] <- adjustcolor(tort.plot.cols2[i],last.params$admix.proportions[i])
}

png(file="tort_admixlocation_spacemix.png",res=200,height=5*200,width=15*200)
	par(mfrow=c(1,3))
	plot(tort.coords,type='n',
			xlab="Long",ylab="Lat",
			main="Tortoise Locations")
		text(tort.coords,labels=tort.names,col=tort.plot.cols2)
	plot(last.params$population.coordinates[1:k,],type='n',
			ylab="'lat'",xlab="'long'",
			main="Tortoise SpaceMix Locations")
		text(last.params$population.coordinates[1:k,],labels=tort.names,col=tort.plot.cols2)
	plot(last.params$population.coordinates[1:k,],type='n',
			xlim=c(-120,-110),ylim=c(32,40),
			ylab="'lat'",xlab="'long'",
			main="Tortoise SpaceMix Locations")
		text(last.params$population.coordinates[1:k,],labels=tort.names,col=tort.plot.cols2)
			points(last.params$population.coordinates[(k+1):(2*k),],
					col= point.tort.ad.cols,pch=20)
	arrows(x0 = last.params$population.coordinates[(k+1):(2*k),1],
			y0 = last.params$population.coordinates[(k+1):(2*k),2],
			x1 = last.params$population.coordinates[1:k,1],
			y1 = last.params$population.coordinates[1:k,2],
			lwd = last.params$admix.proportions,
			col = point.tort.ad.cols,
			length=0.1)
dev.off()






