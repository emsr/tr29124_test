
#include <math.h>


#include "nric.h"


/*
 *    Polynomial coefficient shift.  Given a coefficient array c[0..n-1], this routine
 *    generates a coefficient array d such that  S c_k.y^k = S d_k.x^k  where x and y
 *    are related by y = (x-(b+a)/2)/((b-a)/2).  i.e. the interval -1 < y < +1 is
 *    mapped onto the interval a < x < b.  The array d is returned in c.
 */
void poly_shift_coeff( double a, double b, double *c, int n ) {

    int j, k;
    double factor, con;

    if ( a == b ) nrerror( "Polynomial limits equal in poly_shift_coeff." );
    con = 2.0/(b-a);
    factor = con;
    /*
     *    First, we rescale by the factor con ...
     */
    for ( j = 1; j < n; ++j ) {
        c[j] *= factor;
        factor *= con;
    }
    /*
     *    ... which is then redefined as the desired shift.
     */
    con = (b+a)/2.0;
    /*
     *    We accomplish the shift by synthetic division.
     */
    for ( j = 0; j <= n-2; ++j ) for ( k = n-2; k >= j; --k ) c[k] -= con*c[k+1];
}


