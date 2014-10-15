# compute for each N, M, the chance of picking different alleles
#  and record for each combination of coverages

{ 
    for ( N=1; N<=(NF/2); N++ ) { 
        D[N] = ( $(2*N-1) + $(2*N)  );  # total coverage
        if ( D[N] > 0 && D[N] <= 12 ) {
            P[N] = $(2*N-1) / D[N];  # major allele freq
            if ( D[N] > 1 ) {
                W[N,N,D[N],D[N]] += 1 ;
                PI[N,N,D[N],D[N]] += 2 * P[N] * ( D[N] - $(2*N-1) ) / (D[N]-1) ;  # prob of difference
            }
            for ( M=1; M<N; M++ ) {
                if ( D[M] > 0 && D[M] <= 12 ) {
                    W[N,N,D[N],D[M]] += 1 ;
                    PI[N,M,D[N],D[M]] += ( P[N]*(1-P[M]) + P[M]*(1-P[N]) );
                }
            }
        }
    }
}

END { 
    # outputs mean pi by coverage and total number of sites that fell into that bin
    # and "-1" means NA
    printf "ind1 ind2 " ;
    for ( K=1; K<=12; K++ ){
        for ( L=1; L<=12; L++ ){
            printf "p%i.%i t%i.%i", K, L, K, L;
        } 
    } 
    printf "\n" ;
    for ( N=1; N<=(NF/2); N++ ) { 
        for ( M=1; M<=N; M++ ) { 
            printf "%i %i ", N, M ;
            for ( K=1; K<=12; K++ ){
                for ( L=1; L<=12; L++ ){
                    printf "%f %f ", ( W[N,M,K,L]>0 ? PI[N,M,K,L]/W[N,M,K,L] : -1 ), W[N,M,K,L];
                }
            }
            printf "\n"
        } 
    } 
}
