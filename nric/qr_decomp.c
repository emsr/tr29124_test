

#include <math.h>

#include "nric.h"


void qr_decomp( double **a, int n, double *c, double *d, int *sing ) {

    int i, j, k;
    double scale, sigma, sum, tau;

    *sing = 0.0;
    for ( k = 1; k < n; ++k ) {
        scale = 0.0;
        for ( i = k; i <= n; ++i ) scale = dmax( scale, fabs(a[i][k]) );
        if ( scale == 0.0 ) {
            *sing = 1;
            c[k] = d[k] = 0.0;
        } else {
            for ( i = k; i <= n; ++i ) a[i][k] /= scale;
            for ( sum = 0.0, i = k; i <= n; ++i ) sum += dsqr(a[i][k]);
            sigma = dsign( sqrt(sum), a[k][k] );
            a[k][k] += sigma;
            c[k] = sigma*a[k][k];
            d[k] = -scale*sigma;
            for ( j = k+1; j <= n; ++j ) {
                for ( sum = 0.0, i = k; i <= n; ++i ) sum += a[i][k]*a[i][j];
                tau = sum/c[k];
                for ( i = k; i <= n; ++i ) a[i][j] -= tau*a[i][k];
            }
        }
    }
    d[n] = a[n][n];
    if ( d[n] == 0.0 ) *sing = 1;
}


void qr_backsub( double **a, int n, double *c, double *d, double *b ) {

    int i, j;
    double sum, tau;

    for ( j = 1; j < n; ++j ) {
        for ( sum = 0.0, i = j; i <= n; ++i ) sum += a[i][j]*b[i];
        tau = sum/c[j];
        for ( i = j; i <= n; ++i ) b[i] -= tau*a[i][j];
    }
    r_backsub( a, n, d, b );
}


void r_backsub( double **a, int n, double *d, double *b ) {

    int i, j;
    double sum;

    b[n] /= d[n];
    for ( i = n-1; i >= 1; --i ) {
        for ( sum = 0.0, j = i+1; j <= n; ++j ) sum += a[i][j]*b[j];
        b[i] = ( b[i] - sum )/d[i];
    }
}


void qr_update( double **r, double **qt, int n, double *u, double *v ) {

    int i, j, k;

    for ( k = n; k >= 1; --k ) if ( u[k] ) break;
    if ( k < 1 ) k = 1;
    for ( i = k-1; i >= 1; --i ) {
        jacobi_rotate( r, qt, n, i, u[i], -u[i+1] );
        if ( u[i] == 0.0 ) u[i] = fabs(u[i+1]);
        else if ( fabs(u[i]) > fabs(u[i+1]) ) u[i] = fabs(u[i])*sqrt(1.0 + dsqr(u[i+1]/u[i]));
        else u[i] = fabs(u[i+1])*sqrt(1.0 + dsqr(u[i]/u[i+1]));
    }
    for ( j = 1; j <= n; ++j ) r[1][j] += u[1]*v[j];
    for ( i = 1; i < k; ++i ) jacobi_rotate( r, qt, n, i, r[i][i], -r[i+1][i] );
}


void jacobi_rotate( double **r, double **qt, int n, int i, double a, double b ) {

    int j;
    double c, fact, s, w, y;

    if ( a == 0.0 ) {
        c = 0.0;
        s = dsgn(b);
    } else {
        fact = dpythag( a, b );
        c = a/fact;
        s = b/fact;
    }
    for ( j = i; j <= n; ++j ) {
        y = r[i][j];
        w = r[i+1][j];
        r[i][j] = c*y - s*w;
        r[i+1][j] = s*y + c*w;
    }
    for ( j = 1; j <= n; ++j ) {
        y = qt[i][j];
        w = qt[i+1][j];
        qt[i][j] = c*y - s*w;
        qt[i+1][j] = s*y + c*w;
    }
}


