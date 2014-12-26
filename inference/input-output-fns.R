write.full.hts <- function ( hts, file ) {
    # write out hitting times as a matrix
    colnames( hts ) <- locs
    cat("Writing output to ", file, " .\n")
    write.table( hts, file=file, row.names=FALSE )
}

read.full.hts <- function (file,locs) {
    # file should be in the following format:
    #   "144" "23" "201"
    #    23.0 45.1 120.2
    #    33.3 31.3 110.8
    #    42.0 85.8 124.2
    # and will check the first row corresponds to 'locs'
    hts.locs <- as.numeric( scan(file,nlines=1,what="char") )
    stopifnot( all( hts.locs == locs ) )
    matrix( scan( file, skip=1 ), ncol=length(hts.locs), byrow=TRUE )
}

write.sub.hts <- function ( hts, file ) {
    hts.df <- data.frame(
            row=as.vector(row(hts)), # [upper.tri(hts,diag=TRUE)],
            col=as.vector(col(hts)), # [upper.tri(hts,diag=TRUE)],
            DISTANCE=as.vector(hts)  # [upper.tri(hts,diag=TRUE)]
        )
    cat("Writing out to ", file, " .\n")
    write.table( noisy.df, file=file, row.names=FALSE )
}

read.sub.hts <- function ( file, locs ) {
    # Read in the three-column format of hitting times with row and column labels
    #  TO-DO: should be using tortoise IDs, not row and column numbers like currently.
    ht.df <- read.table( "test_six_layers/256x/six-raster-list-sim-0_00-hts.tsv",header=TRUE,stringsAsFactors=FALSE)
    ht <- matrix( NA, nrow=length(locs), ncol=length(locs) )
    ht[ cbind( ht.df$row, ht.df$col ) ] <- ht.df$DISTANCE
    return( ht )
}
