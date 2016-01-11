

#include <math.h>

#include "nric.h"

void cholesky_decomp( double **a, int n, double *d ) {

    int i, j, k;
    double sum;

    for ( i = 1; i <= n; ++i ) {
        for ( j = 1; j <= n; ++j ) {
            for ( sum = a[i][j], k = i-1; k >= 1; --k ) sum -= a[i][k]*a[j][k];
            if ( i == j ) {
                if ( sum <= 0.0 ) nrerror( "Failure in cholesky_decomp." );
                d[i] = sqrt(sum);
            } else a[j][i] = sum/d[i];
        }
    }
}


void cholesky_backsub( double **a, int n, double *d, double *b, double *x ) {

    int i, k;
    double sum;

    for ( i = 1; i <= n; ++i ) {
        for ( sum = b[i], k = i-1; k >= 1; --k ) sum -= a[i][k]*x[k];
        x[i] = sum/d[i];
    }
    for ( i = n; i >= 1; --i ) {
        for ( sum = x[i], k = i+1; k <= n; ++k ) sum -= a[k][i]*x[k];
        x[i] = sum/d[i];
    }
}

void cholesky_invert( double **a, int n, double *d ) {

    int i, j, k;
    double sum;

    for ( i = 1; i <= n; ++i ) {
        a[i][i] = 1.0/d[i];
        for ( j = i+1; j <= n; ++j ) {
            sum = 0.0;
            for ( k = i; k < j; ++k ) sum -= a[j][k]*a[k][i];
            a[j][i] = sum/d[j];
        }
    }
}



