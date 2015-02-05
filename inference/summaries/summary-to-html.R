#!/usr/bin/Rscript

x <- read.csv("all-collated.csv",stringsAsFactors=FALSE)
x$file <- file.path("..",x$file)
x$fit.png <- paste(gsub("inference-","comparison-",gsub(".RData","-",x$file)),x$summary,".png",sep="")
x$fit.html <- gsub(".RData",".html",x$file)
x$alt_pref <- file.path(dirname(x$file),paste("evaluate-alt_pref_pda-",x$summary,".html",sep=''))
x$alt_1 <- file.path(dirname(x$file),paste("evaluate-alt_1_pda-",x$summary,".html",sep=''))
x$alt_2 <- file.path(dirname(x$file),paste("evaluate-alt_2_pda-",x$summary,".html",sep=''))
x$alt_3 <- file.path(dirname(x$file),paste("evaluate-alt_3_pda-",x$summary,".html",sep=''))
x$alt_4 <- file.path(dirname(x$file),paste("evaluate-alt_4_pda-",x$summary,".html",sep=''))

print.vars <- c("summary","mad","mse","converged","n.refs","file")
xt <- x[,print.vars]
for (xn in c("fit.png","fit.html","alt_pref","alt_1","alt_2","alt_3","alt_4")) {
    xt[xn] <- "none"
    xt[file.exists(x[,xn]),xn] <- paste("<a href='",x[file.exists(x[,xn]),xn],"'>link</a>",sep='')
}
null.ind <- match("null",x$summary)
for (xn in c("summary","fit.png")) {
    xt[null.ind,xn] <- paste("<b>",x[null.ind,xn],"</b>",sep="")
}

xt$mse[xt$mse>1e15] <- Inf
xt$mad[xt$mad>1e15] <- Inf
xt <- xt[,c(setdiff(colnames(xt),"file"),"file")]

require(xtable)
xt.mse <- xtable(xt[order(xt$mse),])
xt.mad <- xtable(xt[order(xt$mad),])
colnames(xt.mse)[2] <- "<a href='all-collated-mad.html'>mad</a>"
colnames(xt.mad)[3] <- "<a href='all-collated-mse.html'>mse</a>"

print.xtable(xt.mse,file="all-collated-mse.html",type='html',sanitize.text.function=identity)
print.xtable(xt.mad,file="all-collated-mad.html",type='html',sanitize.text.function=identity)
