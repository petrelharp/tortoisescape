#!/usr/bin/Rscript

###
# This tortoise falls off the edge of the grid
missing.tort.index <- 57

###
torts <- read.csv("../tort_180_info/1st_180_torts.csv",header=TRUE)
nind <- nrow(torts)

stopifnot( torts$Northing[missing.tort.index] == max(torts$Northing) )


# pairwise distances
tort.dist.table <- read.table("../tort_180_info/1st180_pairwise_distances_sorted_redundancy_removed.txt",header=TRUE)
tort.dists <- numeric(nind^2); dim(tort.dists) <- c(nind,nind)
tort.dists[ cbind( match(tort.dist.table$etort1,torts$EM_Tort_ID), match(tort.dist.table$etort2,torts$EM_Tort_ID) ) ] <- tort.dist.table$DISTANCE
tort.dists <- tort.dists + t(tort.dists)

# info
torts <- torts[-missing.tort.index]
tort.dists <- tort.dists[-missing.tort.index,-missing.tort.index]

save(torts,tort.dists,file="torts-info.RData")
