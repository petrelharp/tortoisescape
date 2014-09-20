torts <- read.csv("../1st_180_torts.csv",header=TRUE)
tortlabs <- gsub("etort.","",torts$EM_Tort_ID)
coverages <- read.csv("../coverage_info.csv")
torts$coverage <- coverages$sequence_yield_megabases[match(torts$EM_Tort_ID,coverages$individual)]

sparsecov <- matrix( scan("alleleCounts500kLoci-covmat.txt"), nrow=180 )
ngscov <- matrix( scan("../tortGen/DS2.covar"), nrow=180)
pmat <- diag(nrow(sparsecov)) - 1/nrow(sparsecov)
proj.sparsecov <- pmat %*% sparsecov %*% pmat

sparsepc <- eigen(proj.sparsecov)
ngspc <- eigen(ngscov)

layout(t(1:2))

# compare covariance matrices
plot(as.vector(sparsecov),as.vector(ngscov)); abline(0,1)
plot(as.vector(sparsecov[lower.tri(ngscov,diag=FALSE)]),as.vector(ngscov[lower.tri(ngscov,diag=FALSE)])); abline(0,1)
# ok, NGStools gets normalized covariance matrix

# compare to normalized covariance matrix
plot(as.vector(proj.sparsecov),as.vector(ngscov)); abline(0,1)
plot(as.vector(proj.sparsecov[lower.tri(ngscov,diag=FALSE)]),as.vector(ngscov[lower.tri(ngscov,diag=FALSE)])); abline(0,1)

# compare PCs
plot(eigen(ngscov)$vectors[,1],eigen(proj.sparsecov)$vectors[,1], type='n'); abline(0,1)
text(eigen(ngscov)$vectors[,1],eigen(proj.sparsecov)$vectors[,1], labels=tortlabs)
plot(eigen(ngscov)$vectors[,2],eigen(proj.sparsecov)$vectors[,2], type='n'); abline(0,1)
text(eigen(ngscov)$vectors[,2],eigen(proj.sparsecov)$vectors[,2], labels=tortlabs)

# residuals correlated to coverage?  yes.
pc1lm <- lm( eigen(ngscov)$vectors[,1] ~ eigen(proj.sparsecov)$vectors[,1] )
pc2lm <- lm( eigen(ngscov)$vectors[,2] ~ eigen(proj.sparsecov)$vectors[,2] )
plot( torts$coverage, resid(pc1lm) )
plot( torts$coverage, resid(pc2lm) )



if (FALSE) {
    tmp <- read.table("alleleCounts500kLoci-subsampled-covmat-awk.txt",header=FALSE)[,3]
    subcov <- matrix(NA, nrow=180, ncol=180)
    subcov[upper.tri(subcov,diag=TRUE)] <- tmp
    subcov[lower.tri(subcov,diag=FALSE)] <- t(subcov)[lower.tri(subcov,diag=FALSE)]

    tmp <- read.table("alleleCounts500kLoci-covmat-awk.txt",header=FALSE)[,3]
    awkcov <- matrix(NA, nrow=180, ncol=180)
    awkcov[upper.tri(awkcov,diag=TRUE)] <- tmp
    awkcov[lower.tri(awkcov,diag=FALSE)] <- t(awkcov)[lower.tri(awkcov,diag=FALSE)]
}
