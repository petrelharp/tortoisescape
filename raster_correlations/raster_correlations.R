################################################################
################################################################
#	calculate correlations between all layers
################################################################
################################################################

if(file.exists("/Volumes/cooplab1/tortoises/100x")){
	setwd("/Volumes/cooplab1/tortoises/100x")
}
require(raster)
rasterOptions(tmpdir=".")
options(error=recover)
#	read in raster files
raster.files <- list.files(pattern=".gri")

#	get a vector of raster names from
#		the names of the raster files
raster.names <- gsub("crop_","",raster.files)
	raster.names <- gsub("resampled_","",raster.names)
	raster.names <- gsub("masked_aggregated_100x_","",raster.names)
	raster.names <- gsub(".gri","",raster.names)

x <- length(raster.files)

#	get all the pairwise comparisons
raster.correlations <- vector("list",length=(x*(x-1)/2))
raster.pairs <- combn(1:x,2)

#	go through each pairwise comparison,
#		read in the rasters for each one,
#		and calculate the cell-wise correlation
#		between them
#	then assign the name of the comparison to that 
#		element of the list
for(i in 1:ncol(raster.pairs)){
	cat(i,"\n")
	raster1 <- raster(raster.files[raster.pairs[1,i]])
	raster2 <- raster(raster.files[raster.pairs[2,i]])
	raster.correlations[[i]] <- layerStats(stack(raster1,raster2),stat="pearson",na.rm=TRUE)
	names(raster.correlations)[i] <- sprintf("%s_X_%s",
													raster.names[raster.pairs[1,i]],
													raster.names[raster.pairs[2,i]])
}

	save(raster.names,raster.correlations,file="raster_correlations.Robj")


################################
#	make a data table
#		to summarize correlations
################################

raster.data.table <- matrix(0,nrow=length(raster.names),ncol=length(raster.names)-1)
row.names(raster.data.table) <- raster.names

for(i in 1:nrow(raster.data.table)){
	cat(i,"\t")
	starts <- paste("^",raster.names[i],"_X",sep="")
	ends <- paste("X_",raster.names[i],"$",sep="")
	tmp.index <- c(which(grepl(starts,names(raster.correlations))),
					which(grepl(ends,names(raster.correlations))))
	tmp.correlations <- lapply(
							lapply(raster.correlations[tmp.index],'[[',1),
						'[',1,2)
	tmp.order <- order(abs(unlist(tmp.correlations)),decreasing=TRUE)
	tmp.corr.values <- round(unlist(tmp.correlations)[tmp.order],3)
	tmp.corr.names <- gsub(starts,"",names(tmp.correlations)[tmp.order])
		tmp.corr.names <- gsub(ends,"",tmp.corr.names)
		tmp.corr.names <- gsub("^_","",tmp.corr.names)
		tmp.corr.names <- gsub("_$","",tmp.corr.names)
	raster.data.table[i,] <- paste(tmp.corr.names,tmp.corr.values,sep=": ")
}


write.csv(raster.data.table,file="full_raster_correlation_table.csv",row.names=TRUE)
write.csv(raster.data.table[,1:5],file="small_raster_correlation_table.csv",row.names=TRUE)

################################
#	form correlations into
#		matrix to visualize
################################

# also, drop some rasters that were won't be using:
#	-summer.precip, as it is apparently a dud
#		all values from raster(summer_precip.tif) are NA
#	-pr_ss2, which has a replacement with much less missing data: pr_ss2_st
#	-bd_ss2, which has a replacement with much less missing data: bd_ss2_st
#	-bedrock_ss2, which has a replacement with much less missing data: bedrock_ss2_st

raster.correlations <- raster.correlations[-grep("summer_precip",names(raster.correlations))]
	starts <- paste("^","pr_ss2","_X",sep="")
	ends <- paste("X_","pr_ss2","$",sep="")
raster.correlations <- raster.correlations[-grep(starts,names(raster.correlations))]
raster.correlations <- raster.correlations[-grep(ends,names(raster.correlations))]
	starts <- paste("^","bd_ss2_30","_X",sep="")
	ends <- paste("X_","bd_ss2_30","$",sep="")
raster.correlations <- raster.correlations[-grep(starts,names(raster.correlations))]
raster.correlations <- raster.correlations[-grep(ends,names(raster.correlations))]
	starts <- paste("^","bedrock_ss2","_X",sep="")
	ends <- paste("X_","bedrock_ss2","$",sep="")
raster.correlations <- raster.correlations[-grep(starts,names(raster.correlations))]
raster.correlations <- raster.correlations[-grep(ends,names(raster.correlations))]
	starts <- paste("^","tree_err_30","_X",sep="")
	ends <- paste("X_","tree_err_30","$",sep="")
raster.correlations <- raster.correlations[-grep(starts,names(raster.correlations))]
raster.correlations <- raster.correlations[-grep(ends,names(raster.correlations))]




raster.names <- raster.names[-grep("summer_precip",raster.names)]
raster.names <- raster.names[-grep(paste("pr_ss2","$",sep=""),raster.names)]
raster.names <- raster.names[-grep(paste("bd_ss2_30","$",sep=""),raster.names)]
raster.names <- raster.names[-grep(paste("bedrock_ss2","$",sep=""),raster.names)]
raster.names <- raster.names[-grep(paste("tree_err_30","$",sep=""),raster.names)]

x <- length(raster.names)
correlation.matrix <- matrix(0,nrow=x,ncol=x)

row.names(correlation.matrix) <- raster.names
colnames(correlation.matrix) <- raster.names

for(i in 1:x){
	cat(i,"\t")
	for(j in 1:x){
		if(i == j){
			correlation.matrix[i,j] <- 1
		} else {
			index <- intersect(
						grep(raster.names[i],names(raster.correlations)),
						grep(raster.names[j],names(raster.correlations)))
			correlation.matrix[i,j] <- raster.correlations[[index]][[1]][1,2]
		}
	}
}

save(correlation.matrix,file="raster_correlation_matrix.Robj")
#then visualize the result
#	first using the raw correlation values,
#	then the abs(corr.values)

#Raw correlation values first

#use heatmap to get the cluster order index
tmp.heatmap <- heatmap(correlation.matrix)

labels <- raster.names[tmp.heatmap$rowInd]
x.ticks <- (1:nrow(correlation.matrix)-1)/(nrow(correlation.matrix)-1)
png(file="raw_correlation_matrix_heatmap.png",res=200,width=10*200,height=10*200)
	par(oma=c(1,1,1,3))
	image(correlation.matrix[tmp.heatmap$rowInd,tmp.heatmap$rowInd],xaxt='n',yaxt='n',main="raster correlations")
		axis(side=1,at=x.ticks,labels=labels,las=2,cex.axis=0.6)
		axis(side=2,at=x.ticks,labels=labels,las=2,cex.axis=0.6)
				leg.x <- c(1:nrow(correlation.matrix))
				leg.y <- sort(as.matrix(correlation.matrix))
				leg.z <- outer(leg.x,leg.y)
				leg.lab.vals <- c(round(min(correlation.matrix),2),round(max(correlation.matrix),2))
				fields::image.plot(leg.x,leg.y,leg.z,add=TRUE,
						nlevel=nrow(correlation.matrix),
						horizontal=FALSE,graphics.reset=TRUE,
						legend.only=TRUE,col=heat.colors(nrow(correlation.matrix)),
						smallplot=c(0.975,0.990,0.25,0.68),				
						axis.args=list(at=c(-82,84),labels=leg.lab.vals))
dev.off()

#Now abs(correlation values)

#use heatmap to get the cluster order index
tmp.heatmap <- heatmap(abs(correlation.matrix))

labels <- raster.names[tmp.heatmap$rowInd]
x.ticks <- (1:nrow(correlation.matrix)-1)/(nrow(correlation.matrix)-1)
png(file="abs_correlation_matrix_heatmap.png",res=200,width=10*200,height=10*200)
	par(oma=c(1,1,1,3))
	image(abs(correlation.matrix)[tmp.heatmap$rowInd,tmp.heatmap$rowInd],xaxt='n',yaxt='n',main="raster correlations")
		axis(side=1,at=x.ticks,labels=labels,las=2,cex.axis=0.6)
		axis(side=2,at=x.ticks,labels=labels,las=2,cex.axis=0.6)
				leg.x <- c(1:nrow(correlation.matrix))
				leg.y <- sort(as.matrix(correlation.matrix))
				leg.z <- outer(leg.x,leg.y)
				leg.lab.vals <- c(round(min(abs(correlation.matrix)),2),round(max(abs(correlation.matrix)),2))
				fields::image.plot(leg.x,leg.y,leg.z,add=TRUE,
						nlevel=nrow(correlation.matrix),
						horizontal=FALSE,graphics.reset=TRUE,
						legend.only=TRUE,col=heat.colors(nrow(correlation.matrix)),
						smallplot=c(0.975,0.990,0.25,0.68),				
						axis.args=list(at=c(-81,83),labels=leg.lab.vals))
dev.off()


#	and, finally, some PCs w/ the correlation matrix
abs.correlation.matrix <- abs(correlation.matrix)
corr.eigen.obj <- eigen(abs.correlation.matrix)
PC.loadings <- corr.eigen.obj$values/sum(corr.eigen.obj$values)

for(i in 1:10){
	pdf(file=paste("PC_",i,".pdf",sep=""),width=7,height=7)
		plot(corr.eigen.obj$vectors[tmp.heatmap$rowInd,i],type='n',
				main=paste("PC_",i,"\n","Var. explained:",round(PC.loadings[i],3)))
			text(corr.eigen.obj$vectors[tmp.heatmap$rowInd,i],
					labels=row.names(correlation.matrix)[tmp.heatmap$rowInd],cex=0.5)
	dev.off()
}

if(FALSE){
	plot(corr.eigen.obj$vectors[,1],type='n')
		text(corr.eigen.obj$vectors[,1],labels=row.names(correlation.matrix),cex=0.5)
	plot(corr.eigen.obj$vectors[,1],corr.eigen.obj$vectors[,2],type='n')
		text(corr.eigen.obj$vectors[,1],corr.eigen.obj$vectors[,2],labels=row.names(correlation.matrix),cex=0.5)
}




