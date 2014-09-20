
sparsecov <- matrix( scan("alleleCounts_1millionloci-subsampled-covmat.txt"), nrow=180 )
ngscov <- matrix( scan("../tortGen/DS2.covar"), nrow=180)

layout(t(1:2))
plot(as.vector(sparsecov),as.vector(ngscov)); abline(0,1)
plot(as.vector(sparsecov[lower.tri(ngscov,diag=FALSE)]),as.vector(rcov[lower.tri(subcov,diag=FALSE)])); abline(0,1)


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
