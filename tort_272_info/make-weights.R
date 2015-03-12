#!/usr/bin/Rscript

require(raster)

dups <- list(
        c("etort-156","etort-296"),
        c("etort-143","etort-297")
    )

# geographic distance on which weights decay (meters)
weight.scale <- 25e3

sample.info <- read.csv("sample_metadata.csv",header=TRUE,stringsAsFactors=FALSE)
geodist.tab <- read.csv("geog_distance.csv",header=TRUE,stringsAsFactors=FALSE)
geodist <- matrix(NA,nrow=nrow(sample.info),ncol=nrow(sample.info))
dimnames(geodist) <- list( sample.info$EM_Tort_ID, sample.info$EM_Tort_ID )
geodist.inds <- cbind( match(geodist.tab[,1],rownames(geodist)), match(geodist.tab[,2],colnames(geodist)) )
usethese <- apply( !is.na(geodist.inds), 1, all )
geodist[ geodist.inds[usethese,] ] <- geodist.tab[usethese,3]
geodist[is.na(geodist)] <- t(geodist)[is.na(geodist)]

weights <- 1/rowSums( exp(-geodist/weight.scale) )
for (dup in dups) {
    dup <- dup[ dup %in% names(weights) ]
    cat("Downweighting ", paste(dup,collapse=" and "), ".\n")
    weights[ dup ] <- mean(weights[dup])/length(dup)
}

weight.tab <- data.frame(etort=names(weights),weight=weights)

write.csv(weight.tab,file="weights.csv",row.names=FALSE)

# also, one that's just downweighted by duplicates
weights.nodup <- weights
weights.nodup[] <- 1
for (dup in dups) {
    dup <- dup[ dup %in% names(weights.nodup) ]
    weights.nodup[ dup ] <- mean(weights.nodup[dup])/length(dup)
}
weight.nodup.tab <- data.frame(etort=names(weights.nodup),weight=weights.nodup)

write.csv(weight.nodup.tab,file="weights-nodups.csv",row.names=FALSE)
