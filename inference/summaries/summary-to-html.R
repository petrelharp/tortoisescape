#!/usr/bin/Rscript

x <- read.csv("all-collated.csv",stringsAsFactors=FALSE)
x$file <- file.path("..",x$file)
x$fit.all.png <- gsub("inference-","comparison-",gsub(".RData","-all.png",x$file))
x$fit.habitat.png <- gsub("inference-","comparison-",gsub(".RData","-habitat-only.png",x$file))
x$eval.all <- file.path(dirname(x$file),"evaluate-alt_pref_pda-all.html")
x$eval.habitat <- file.path(dirname(x$file),"evaluate-alt_pref_pda-habitat-only.html")

print.vars <- c("summary","mad","mse","converged","n.refs","file")
xt <- ( cbind( x[,print.vars], 
                fit.all.png="none",
                fit.habitat.png="none",
                eval.all="none",
                eval.habitat="none"
            , stringsAsFactors=FALSE) )
for (xn in c("fit.all.png","fit.habitat.png","eval.all","eval.habitat")) {
    xt[file.exists(x[,xn]),xn] <- paste("<a href='",x[file.exists(x[,xn]),xn],"'>link</a>",sep='')
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
