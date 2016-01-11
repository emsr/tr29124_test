
#include <math.h>

#include "nric.h"


/*
 *    Returns the derivative of a function funk at a point x by Ridder's method of polynomial extrapolation.
 *    The valie h is input as an estimated stepsix=ze;  it should not be small but rather it should be an
 *    interval over which funk changes substantially.  An estmate of the error in the derivative is returned
 *    as err.
 */
double diff_ridder( double (*funk)(double), double x, double h, double *err ) {

    int i, j;
    double errt, fac, hh, **a, ans;

    const int NTAB = 10;
    const double CON = 1.4, CON2 = 1.4*1.4, BIG = 1.0e30, SAFE = 2.0;

    if ( h <= 0.0 ) nrerror( "Stepsize must be posative in diff_ridder." );

    a = dmatrix( 1, NTAB, 1, NTAB );
    hh = h;
    a[1][1] = ((*funk)(x+hh) - (*funk)(x-hh))/(2.0*hh);
    *err = BIG;
    /*
     *    Successive columns of the Neville tableau will go to smaller stepsizes
     *    and to higher orders of extrapolation.
     */
    for ( j = 2; j <= NTAB; ++j ) {
        /*
         *    try a new, smaller stepsize.
         */
        hh /= CON;
        a[1][j] = ((*funk)(x+hh) - (*funk)(x-hh))/(2.0*hh);
        fac = CON2;
        for ( i = 2; i <= j; ++i ) {
            /*
             *    Compute extrapolations of various orders, requiring
             *    no new function evaluations.
             */
            a[i][j] = (fac*a[i-1][j] - a[i-1][j-1])/(fac - 1.0);
            fac *= CON2;
            /*
             *    The error strategy is to compare the each new extrapolation to one order lower,
             *    both at the present stepsize and to the previous one.
             */
            errt = dmax(fabs(a[i][j] - a[i-1][j]), fabs(a[i][j] - a[i-1][j-1]) );
            if ( errt <= *err ) {
                *err = errt;
                ans = a[i][j];
            }
        }
        /*
         *    If higher order is worse by a significant factor SAFE then quit early.
         */
        if ( fabs(a[j][j]) - fabs(a[j-1][j-1]) >= SAFE*(*err) ) break;
    }
    free_dmatrix( a, 1, NTAB, 1, NTAB );
    return ans;
}


