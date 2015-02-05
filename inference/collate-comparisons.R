#!/usr/bin/Rscript

usage <- "Collate results produced by comparison-results.R, stored in the .RData files passed in as arguments.
Usage:
    Rscript collate-comparisons.R (outfile) ( file names ) 
"

argvec <- if (interactive()) { scan(what='char') } else { commandArgs(TRUE) }
if (length(argvec) < 2) { stop(usage) }

outfile <- argvec[1]
infiles <- argvec[-1]
readable <- file.exists(infiles)
for (k in which(!readable)) { warning(infiles[k], " does not exist.\n") }

source("resistance-fns.R")
require(raster)

gmat <- function (geodist.tab,pimat) {
    geodist <- pimat
    geodist[] <- NA
    geodist.inds <- cbind( match(geodist.tab[,1],rownames(geodist)), match(geodist.tab[,2],colnames(geodist)) )
    usethese <- apply( !is.na(geodist.inds), 1, all )
    geodist[ geodist.inds[usethese,] ] <- geodist.tab[usethese,3]
    geodist[is.na(geodist)] <- t(geodist)[is.na(geodist)]
    geodist
}

# null model fit
null.config.file <- "summaries/all/config.json"
null.config <- read.json.config(null.config.file)
null.env <- new.env()
load(file.path(dirname(null.config.file),null.config$setup_files),envir=null.env)
assign("geodist.tab", read.csv( file.path(dirname(null.config.file),dirname(null.config$sample_locs),"geog_distance.csv"), header=TRUE, stringsAsFactors=FALSE ), null.env )
assign("pcs", read.csv(file.path(dirname(null.config.file),dirname(null.config$divergence_file),"pcs.csv"),header=TRUE), null.env )
assign("geodist", with(null.env, { gmat(geodist.tab,pimat) } ), null.env )
null.results <- with( null.env, {
            nearby.weights <- 1 / rowSums( geodist < 25e3 )
            pairwise.weights <- outer(nearby.weights,nearby.weights,"*")
            omit.comparisons <- ( pcs$PC1[match(rownames(pimat),pcs$etort)][row(pimat)] * pcs$PC1[match(colnames(pimat),pcs$etort)][col(pimat)] < 0 )
            dup.inds <- match( c( "etort-296", "etort-297" ), rownames(pimat) )
            omit.comparisons <- ( omit.comparisons | (row(pimat) %in% dup.inds) | (col(pimat) %in% dup.inds) )
            # and omit self comparisons and ONLY UPPER TRIANGLE
            omit.comparisons <- ( omit.comparisons | (row(pimat) < col(pimat))  )
            resids <- resid( lm( pimat[!omit.comparisons] ~ geodist[!omit.comparisons] ) )
            # weighted median abs( residual )
            w.mad <- weighted.median( abs(resids), pairwise.weights[!omit.comparisons] )
            w.mse <- sqrt( weighted.mean( resids^2, pairwise.weights[!omit.comparisons], na.rm=TRUE ) )
            list( 
                summary="null",
                mad=w.mad, 
                mse=w.mse, 
                converged=NA,
                n.refs=NA,
                file=NA
            )
        } )


results <- c( list(null.results), lapply( infiles[readable], function (infile) {
        load(infile)
        pcs <- read.csv(file.path(dirname(config.file),dirname(config$divergence_file),"pcs.csv"),header=TRUE)
        omit.comparisons <- ( pcs$PC1[match(rownames(pimat),pcs$etort)][row(pimat)] * pcs$PC1[match(colnames(pimat),pcs$etort)][col(pimat)] < 0 )
        pc.cols <- adjustcolor( ifelse( pcs$PC1[match(rownames(pimat),pcs$etort)][row(pimat)] < 0, "purple", "blue" ), 0.5 )
        # remove duplicates: these are (etort-156 / etort-296 ) and (etort-143 / etort-297)
        dup.inds <- match( c( "etort-296", "etort-297" ), rownames(pimat) )
        omit.comparisons <- ( omit.comparisons | (row(pimat) %in% dup.inds) | (col(pimat) %in% dup.inds) )
        # and omit self comparisons
        omit.comparisons <- ( omit.comparisons | (row(pimat) == col(pimat))  )

        if (is.numeric(hts)) {
            fitted <- paramvec(local.config)[1] + (hts+t(hts))/2 
            resids <- (fitted - pimat)
            resids[omit.comparisons] <- NA
            fitted[omit.comparisons] <- NA

            # weight residuals by 1 / number of other samples within 25km
            geodist <- gmat(geodist.tab,pimat)
            nearby.weights <- 1 / rowSums( geodist < 25e3 )
            pairwise.weights <- outer(nearby.weights,nearby.weights,"*")

            ut <- upper.tri(pimat,diag=FALSE)
            # weighted median abs( residual )
            w.mad <- weighted.median( abs(resids)[ut], pairwise.weights[ut] )
            w.mse <- sqrt( weighted.mean( resids[ut]^2, pairwise.weights[ut], na.rm=TRUE ) )
        } else {
            w.mad <- w.mse <- NA
        }

        return( list(
                summary=basename(dirname(config.file)), 
                mad=w.mad, 
                mse=w.mse, 
                converged=trust.optim.results$converged, 
                n.refs=length(local.config$reference_inds), 
                file=infile 
            ) )
    } ) )

cat("Saving results to ", outfile, "\n")

save(results, file=outfile)

csvfile <- gsub("RData$","csv",outfile)
cat(" and writing to ", csvfile, "\n")

results.tab <- do.call(rbind,lapply(results,as.data.frame))
results.tab <- results.tab[order(results.tab[,"mad"]),]

write.csv( results.tab, file=csvfile, row.names=FALSE )
