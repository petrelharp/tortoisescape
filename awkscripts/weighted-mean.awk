# compute means of odd rows
# weighted by sum of odd and even rows

NR > 1 { 
    for ( N=1; N<=(NF/2); N++ ) { 
        A[N]+= ( $(2*N-1) )*( $(2*N-1) + $(2*N)  )
        W[N]+= $(2*N)
        for ( M=N; M<=NF; M++ ) {
            C[N,M] += $N * $M
        }
    }
}

END { for ( N=1; N<=NF/2; N++ ) print A[N]/W[N] }
