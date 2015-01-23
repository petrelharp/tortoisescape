#!/usr/bin/Rscript

usage <- "Write PCA info from covariance matrix.
Usage:
   Rscript write-PCs.R (covariance table) (output file)
input file is in the format:
 etort1 etort2 covariance
  ...
"

argvec <- if (interactive()) { scan(what="char") } else { commandArgs(TRUE) }

if (length(argvec)<2) { stop(usage) }
covmat.file <- argvec[1]
outfile <- argvec[2]

## PCA
cov.table <- read.csv(covmat.file,header=TRUE,stringsAsFactors=FALSE)
inds <- (unique(c(cov.table[,1],cov.table[,2])))
covmat <- matrix(NA,nrow=length(inds),ncol=length(inds))
dimnames(covmat) <- list(inds,inds)
covmat[ cbind( match(cov.table[,1],inds), match(cov.table[,2],inds) ) ] <- cov.table[,3]
covmat[is.na(covmat)] <- t(covmat)[is.na(covmat)]

pmat <- diag(length(inds)) - 1/length(inds)
eig.covmat <- eigen(pmat %*% covmat %*% pmat)

# proportion of variance explained
# eig.covmat$values / sum(eig.covmat$values)
cat( "Proportion of variance:\n" )
cat( paste(formatC( eig.covmat$values / sum(eig.covmat$values)*100, digits=3 )[1:8],"%",sep=''), "\n" )

pcs <- eig.covmat$vectors[,1:10]
colnames(pcs) <- paste("PC",1:ncol(pcs),sep='')
pcs <- cbind( data.frame( etort=inds ), pcs )

write.csv( pcs, file=outfile, row.names=FALSE )
