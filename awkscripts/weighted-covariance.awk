# compute covariances of ratio of odd rows to sum of odd and even rows
# weighted by geometric mean of corresponding even rows

NR == 1 { for ( N=1; N<=NF; N++ ) H[N] = $N }

NR > 1 { 
    for ( N=1; N<=(NF/2); N++ ) { 
        D[N] = ( $(2*N-1) + $(2*N)  );  # total coverage
        P[N] = $(2*N-1) / D[N];  # major allele freq
        for ( M=1; M<=N; M++ ) {
            WW = sqrt(D[N]*D[M]);
            W[N,M] += WW;
            C[N,M] += ( P[N] )*( P[M] )*WW;
            Z[N,M] += ( P[N] )*WW;
            if (N!=M) { Z[M,N] += ( P[M] )*WW;
        }
    }
}

END { for ( N=1; N<=(NF/2); N++ ) { for ( M=1; M<=N; M++ ) { print H[N],H[M],C[N,M]/W[N,M]-(Z[N,M]/W[N,M])*(Z[M,N]/W[N,M]); } } }
