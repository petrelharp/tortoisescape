{ 
    for (N=1; N<=(NF/3); N++) {
        R[N]=$(3*N-2);
        H[N]=$(3*N-1);
        A[N]=$(3*N);
        PI[N,N] += H[N];
        for (M=1; M<N; M++) {
            # pi = aR*bA + 0.5*(aR*bH + aH*bR + aH*bA + aA*bH + aH*bH) + aA*bR
            PI[N,M] += R[N]*A[M] + 0.5 * ( R[N]*H[M] + H[N]*R[M] + A[N]*H[M] + H[N]*A[M] + H[N]*H[M] ) + A[N]*R[M] ;
        }
    }
}

END {
    for (N=1; N<=(NF/3); N++) {
        for (M=1; M<=N; M++) {
            print PI[N,M]/NR;
        }
    }
}
