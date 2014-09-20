# compute covariances of nonzeroness of odd rows
# at those rows where the sum of odd and even rows are nonzero

NR == 1 { srand(); for ( N=1; N<=NF; N++ ) H[N] = $N }

NR > 1 { 
    for ( N=1; N<=(NF/2); N++ ) { 
        D[N] = ( $(2*N-1) + $(2*N)  );  # total coverage
        if ( D[N] > 0 ) {
            A[N] = rand() < $(2*N-1)/D[N] ? 1 : 0 ;  # pick a random allele
            for ( M=1; M<=N; M++ ) {
                # print NR,":",N,M,$(2*N),$(2*M);
                if ( $(2*M) > 0 ) {
                    W[N,M]++;  # number of sites where both have coverage
                    C[N,M] += A[N]*A[M];
                    Z[N,M] += A[N];
                    if (N!=M) { Z[M,N] += A[M]; }
                    # print "W:",W[N,M],"C:",C[N,M],"Z:",Z[M,N],Z[N,M];
                }
            }
        }
    }
}

END { for ( N=1; N<=(NF/2); N++ ) { for ( M=1; M<=N; M++ ) { print H[N],H[M],C[N,M]/W[N,M]-(Z[N,M]/W[N,M])*(Z[M,N]/W[N,M]); } } }
