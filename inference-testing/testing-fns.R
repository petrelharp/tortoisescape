##

plot.ht <- function (ht,dims=c(sqrt(length(ht)),sqrt(length(ht)))) {
    dim(ht) <- dims; image(ht)
}

plot.hts <- function (hts,dims=c(sqrt(nrow(hts)),sqrt(nrow(hts)))) {
    for (k in 1:ncol(hts)) {
        plot.ht(hts[,k],dims=dims)
        readline("next?")
    }
}

