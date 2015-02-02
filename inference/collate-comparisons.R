#!/usr/bin/Rscript

usage <- "Collate results stored in the .RData files passed in as arguments.
Usage:
    Rscript collate-comparisons.R ( file names )
"

argvec <- if (interactive()) { scan(what='char') } else { commandArgs(TRUE) }
if (length(argvec)==0) { stop(usage) }

infiles <- argvec
readable <- file.exists(argvec)
for (k in which(!readable)) { warning(infiles[k], " does not exist.\n") }

source("resistance-fns.R")
require(raster)

results <- lapply( infiles[readable], function (infile) {
        load(infile)
        pcs <- read.csv(file.path(dirname(config.file),dirname(config$divergence_file),"pcs.csv"),header=TRUE)
        omit.comparisons <- ( pcs$PC1[match(rownames(pimat),pcs$etort)][row(pimat)] * pcs$PC1[match(colnames(pimat),pcs$etort)][col(pimat)] < 0 )
        pc.cols <- adjustcolor( ifelse( pcs$PC1[match(rownames(pimat),pcs$etort)][row(pimat)] < 0, "purple", "blue" ), 0.5 )
        # remove duplicates: these are (etort-156 / etort-296 ) and (etort-143 / etort-297)
        dup.inds <- match( c( "etort-296", "etort-297" ), rownames(pimat) )
        omit.comparisons <- ( omit.comparisons | (row(pimat) %in% dup.inds) | (col(pimat) %in% dup.inds) )
        # and omit self comparisons
        omit.comparisons <- ( omit.comparisons | (row(pimat) == col(pimat))  )

        fitted <- paramvec(local.config)[1] + (hts+t(hts))/2 
        resids <- (fitted - pimat)
        resids[omit.comparisons] <- NA
        fitted[omit.comparisons] <- NA

        # weight residuals by 1 / number of other samples within 25km
        geodist <- pimat
        geodist[] <- NA
        geodist.tab <- read.csv( file.path(dirname(config.file),dirname(config$sample_locs),"geog_distance.csv"), header=TRUE, stringsAsFactors=FALSE )
        geodist.inds <- cbind( match(geodist.tab[,1],rownames(geodist)), match(geodist.tab[,2],colnames(geodist)) )
        usethese <- apply( !is.na(geodist.inds), 1, all )
        geodist[ geodist.inds[usethese,] ] <- geodist.tab[usethese,3]
        geodist[is.na(geodist)] <- t(geodist)[is.na(geodist)]
        nearby.weights <- 1 / rowSums( geodist < 25e3 )
        pairwise.weights <- outer(nearby.weights,nearby.weights,"*")

        med.resids <- apply(resids,1,weighted.median,w=nearby.weights)

        ut <- upper.tri(pimat,diag=FALSE)
        # weighted median abs( residual )
        w.mad <- weighted.median( abs(resids)[ut], pairwise.weights[ut] )
        w.mse <- sqrt( weighted.mean( resids[ut]^2, pairwise.weights[ut], na.rm=TRUE ) )

        return( 
                summary=basename(dirname(config.file)), 
                mad=w.mad, 
                mse=w.mse, 
                converged=trust.optim.results$converged, 
                n.refs=length(local.config$reference_inds), 
                file=infile 
            )

    } )
