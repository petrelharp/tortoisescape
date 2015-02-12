# compute covariances of ratio of odd rows to sum of odd and even rows
# weighted by product of sums of odd and even rows

NR == 1 { for ( N=1; N<=NF; N++ ) H[N] = $N }  # read in column labels (should be one label for every two rows)

NR > 1 { 
    for ( N=1; N<=(NF/2); N++ ) { 
        D[N] = ( $(2*N-1) + $(2*N)  );  # total coverage
        if ( D[N] > 0 && D[N] <= 12 ) {
            P[N] = $(2*N-1);  # major allele counts
            # diagonal
            if ( D[N] > 1 ) {
                W[N,N] += (D[N]*(D[N]-1));   # weights
                C[N,N] += ( P[N] )*( P[N] )*(D[N]-1)/( D[N] );  # weighted sum of products
                Z[N,N] += ( P[N] )*(D[N]-1);  # weighted sum of first one (note cancelation of D[N])
            }
            for ( M=1; M<N; M++ ) {
                if ( D[M] > 0 && D[M] <= 12 ) {
                    W[N,M] += (D[N]*D[M]);   # weights
                    C[N,M] += ( P[N] )*( P[M] );  # weighted sum of products
                    Z[N,M] += ( P[M] )*D[N];  # weighted sum of first one (note cancelation of D[M])
                    Z[M,N] += ( P[N] )*D[M];  # weighted sum of second one (IN TRANSPOSE)
                }
            }
        }
    }
}

END { 
    for ( N=1; N<=(NF/2); N++ ) { 
        for ( M=1; M<=N; M++ ) { 
            print H[N], H[M], (C[N,M]/W[N,M]) - (Z[N,M]/W[N,M])*(Z[M,N]/W[N,M]); 
        } 
    } 
}
