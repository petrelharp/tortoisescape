# compute means of odd rows
# weighted by even rows

NR == 1 { for ( N=1; N<=NF; N++ ) H[N] = $N }

NR > 1 { 
    for ( N=1; N<=(NF/2); N++ ) { 
        for ( M=1; M<=N; M++ ) {
            W[N,M] += sqrt($(2*N)*$(2*M));
            C[N,M] += ( $(2*N-1) )*( $(2*M-1) )*sqrt($(2*N)*$(2*M));
            Z[N,M] += ( $(2*N-1) )*sqrt($(2*N)*$(2*M));
            if (N!=M) { Z[M,N] += ( $(2*M-1) )*sqrt($(2*N)*$(2*M)); }
        }
    }
}

END { for ( N=1; N<=(NF/2); N++ ) { for ( M=1; M<=N; M++ ) { print H[N],H[M],C[N,M]/W[N,M]-(Z[N,M]/W[N,M])*(Z[M,N]/W[N,M]); } } }
