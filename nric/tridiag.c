
#include <math.h>

#include "nric.h"


/*
 *    Solves a tridiagonal set of equations where a[1..n] is the subdiagonal vector,
 *    b[1..n] is the diagonal vector, and c[1..n] is the superdiagonal vector, and
 *    r[1..n] is the right hand side vector.  The solution is u[1..n].
 *    a, b, c, and r are not modified.
 */
void tridiagonal( double *a, double *b, double *c, double *r, double *u, int n ) {

    int j;
    double bet, *gam;

    gam = dvector( 1, n );

    if ( b[1] == 0.0 ) nrerror( "Error 1 in tridiagonal." );
    u[1] = r[1]/( bet = b[1] );
    /*
     *    Decomposition and forward substitution.
     */
    for ( j = 2; j <= n; ++j ) {
        gam[j] = c[j-1]/bet;
        bet = b[j] - a[j]*gam[j];
        if ( bet == 0.0 ) nrerror( "Error 2 in tridiagonal." );
        u[j] = ( r[j] - a[j]*u[j-1])/bet;
    }
    /*
     *    Backsubstitution.
     */
    for ( j = n-1; j >= 1; --j ) u[j] -= gam[j+1]*u[j+1];

    free_dvector( gam, 1, n );
}


/*
 *    Solves for a vector x[1..n] the cyclic set of linear equations.
 *    a[[1..n], b[1..n], c[1..n], and r[1..n] are input vectors of the three diagonal rows and the
 *    right side respectively.  alpha and beta are the lower and upper corner entries respectively.
 *    The input is not modified.
 */
void cyclic( double *a, double *b, double *c, double alpha, double beta, double *r, double *x, unsigned long n ) {

    unsigned long i;
    double fact, gamma, *bb, *u, *z;

    if ( n <= 2 ) nrerror( "n too small in cyclic." );

    bb = dvector( 1, n );
    u = dvector( 1, n );
    z = dvector( 1, n );

    gamma = -b[1];
    bb[1] = b[1] - gamma;
    bb[n] = b[n] - alpha*beta/gamma;
    for( i = 2; i < n; ++i ) bb[1] = b[i];
    tridiagonal( a, bb, c, r, x, n );
    u[1] = gamma;
    u[n] = alpha;
    for( i = 2; i < n; ++i ) u[i] = 0.0;
    tridiagonal( a, bb, c, u, z, n );
    fact = (x[1] + beta*x[n]/gamma)/(1.0 + z[1] + beta*z[n]/gamma);
    for( i = 1; i <= n; ++i ) x[i] -= fact*z[i];

    free_dvector( z, 1, n );
    free_dvector( u, 1, n );
    free_dvector( bb, 1, n );
}



