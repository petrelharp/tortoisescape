# compute means of columns in input

NR > 1 { for ( N=1; N<=NF; N++ ) A[N]+=$N }

END { for ( N=1; N<=NF; N++ ) print A[N]/(NR-1) }
