
#include <math.h>

#include "nric.h"


/*
 *    Belch out the jacobi polynomials.
 *    Why not!!!
 */

double poly_jacobi( int n, double alpha, double beta, double x ) {

    int j;
    double c, d, e, f, pm1, pp1, p;

    pm1 = 1.0;
    if ( n == 0 ) return pm1;

    p = 0.5*(alpha - beta + (2.0 + alpha + beta)*x);
    if ( n == 1 ) return p;

    for ( j = 1; j < n; ++j ) {

        c = 2.0*(j + 1)*(j + alpha + beta + 1)*(2*j + alpha + beta);
        d = (2*j + alpha + beta+1)*(alpha*alpha - beta*beta);
        e = (2*j + alpha + beta)*(2*j + alpha + beta + 1)*(2*j + alpha + beta + 2);
        f = 2.0*(j + alpha)*(j + beta)*(2*j + alpha + beta + 2);

        if ( c == 0.0 ) nrerror( "Error in poly_jacobi." );
        pp1 = ((d + e*x)*p - f*pm1)/c;
        pm1 = p;
        p = pp1;
    }
    return p;
}


