
#include <math.h>


#include "nric.h"


void poly_interp( double *xa, double *ya, int n, double x, double *y, double *dy ) {

    int i, m, ns = 1;
    double den, dif, dift, ho, hp, w;
    double *c, *d;

    c = dvector( 1, n );
    d = dvector( 1, n );

    dif = fabs(x - xa[1]);

    /*
     *    Here we find the index of the closest table entry,
     */
    for ( i = 1; i <= n; ++i ) {
        if ( ( dift = fabs(x - xa[i]) ) < dif ) {
            ns = i;
            dif = dift;
        }
        /*
         *    And initialize the tableau of c's and d's.
         */
        d[i] = c[i] = ya[i];
    }

    /*
     *    This is the initial approximation to y.
     */
    *y = ya[ns--];

    /*
     *    For each column of the tableau, loop over the current c's and d's and update them.
     */
    for ( m = 1; m < n; ++m ) {
        for ( i = 1; i <= n-m; ++i ) {
            ho = xa[i] - x;
            hp = xa[i+m] - x;
            w = c[i+1] - d[i];
            /*
             *    This error can occur if two input xa's are identical within roundoff error.
             */
            if ( ( den = ho - hp ) == 0.0 ) nrerror( "Error in routine poly_interp." );
            den = w/den;
            /*
             *    Here, the c's and d's are updated.
             */
            d[i] = hp*den;
            c[i] = ho*den;
        }
        *y += ( *dy = ( 2*ns < (n - m) ? c[ns + 1] : d[ns--] ) );
    }
    free_dvector( d, 1, n );
    free_dvector( c, 1, n );
    return;
}



/*
 *    Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns
 *    a value of y and an accuracy estimate dy.  The value returned is that of the 
 *    diagonal rational function, evaluated at x, which passes through the
 *    n points (xa[i], ya[i]), i = 1..n.
 */
void rat_interp( double *xa, double *ya, int n, double x, double *y, double *dy ) {

    int m, i, ns = 1;
    double w, t, hh, h, dd, *c, *d;
    const double TINY = 1.0e-25;

    c = dvector( 1, n );
    d = dvector( 1, n );

    hh = fabs( x - xa[1] );
    for ( i = 1; i <= n; ++i ) {
        h = fabs( x - xa[i] );
        if ( h == 0.0 ) {
            *y = ya[i];
            *dy = 0.0;
            free_dvector( d, 1, n );
            free_dvector( c, 1, n );
            return;
        } else if ( h < hh ) {
            ns = i;
            hh = h;
        }
        c[i] = ya[i];
        /*
         *    The tiny part is needed to prevent a rare zero over zero condition.
         */
        d[i] = ya[i] + TINY;
    }
    *y = ya[ns--];
    for ( m = 1; m < n; ++m ) {
        for ( i = 1; i <= n-m; ++i ) {
            w = c[i+1] - d[i];
            h = xa[i+m] - x;
            t = ( xa[i] - x )*d[i]/h;
            dd = t - c[i+1];
            /*
             *    This error condition indicated that the interpolating function has a pole 
             *    at the requested value of x.
             */
            if ( dd == 0.0 ) nrerror( "Error in routine rat_interp." );
            dd = w/dd;
            d[i] = c[i+1]*dd;
            c[i] = t*dd;
        }
        *y += ( *dy = ( 2*ns < (n-m) ? c[ns+1] : d[ns--] ) );
    }
    free_dvector( d, 1, n );
    free_dvector( c, 1, n );
    return;
}



/**************************************************************************************************************************

                spline()

    Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e. y[i] = f(x[i]), with monatonically increasing
    values of x[i], and given the values yp1, and ypn for the first derivatives of the function at points 1 and n 
    respectively, this routine returns the an array ypp[1..n] that contains the second derivatives of the interpolating 
    function at the tabulated points x[i].  If yp1 and/or ypn are equal to 1.0e30 or larger, the routine is signalled to 
    set the corresponding boundary condition for a natural spline, with zero second derivative at that boundary.

**************************************************************************************************************************/

void  spline( double *x, double *y, int n, double yp1, double ypn, double *ypp ) {

    int i, j;
    double p, qn, sigma, un, *u;

    u = dvector( 1, n - 1 );

    /*
     *    The lower boundary condition is set to either the natural one or
     *    to match a specified first derivative.
     */
    if ( yp1 >= 1.0e30 ) ypp[1] = u[1] = 0.0;
    else {
        ypp[1] = -0.5;
        u[1] = (3.0/(x[2] - x[1]))*((y[2] - y[1])/(x[2] - x[1]) - yp1);
    }

    /*
     *    Decomposition of the tridiagonal system.
     */
    for ( i = 2; i <= n - 1; ++i ) {
        sigma = (x[i] - x[i-1])/(x[i+1] - x[i-1]);
        p = sigma*ypp[i-1] + 2.0;
        ypp[i] = (sigma - 1.0)/p;
        u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1]);
        u[i] = (6.0*u[i]/(x[i+1] - x[i-1]) - sigma*u[i-1])/p;
    }
    /*
     *    The upper boundary condition is set to either the natural one or
     *    to match a specified first derivative.
     */
    if ( ypn >= 1.0e30 ) qn = un = 0.0;
    else {
        qn = 0.5;
        un = (3.0/(x[n] - x[n-1]))*(ypn - (y[n] - y[n-1])/(x[n] - x[n-1]));
    }

    /*
     *    This is the backsubstitution loop of the tridiagonal algorithm.
     */
    ypp[n] = ( un - qn*u[n-1] )/( qn*ypp[n-1] + 1.0 );
    for ( j = n - 1;j >= 1; --j ) ypp[j] = ypp[j]*ypp[j+1] + u[j];

    free_dvector( u, 1, n - 1 );
}




/**************************************************************************************************************************

                spline_interp

        Given the arrays xa[1..na], ya[1..na] which tabulate a function (with the x[i]'s in order), and given the array 
        yapp[1..na] which is the output from spline() above, and given a value of x, this routine returns a cubic-spline 
        interpolated value y.

**************************************************************************************************************************/


double spline_interp( double *xa, double *ya, double *yapp, int na, double x ) {

    int k, klo, khi;
    double h, b, a;

    /*
     *    Store the last values of klo and khi to see if they still work on subsequent calls.
     *    Also store the last number of points to watch out for new spline data.
     *    If the old points don't work reset them to the endpoints.
     *    In either case, refine klo and khi with bisection.
     */
    static int naold = 0;
    static int kloold;
    static int khiold;

    /*
     *    If the number of points has changed then spline was rerun or a new set of points
     *    is being provided.  To prevent the arrays being indexed to an illegal value and
     *    to reflect the arrival of new data reset the values of kloold and khiold and naold.
     */
    if ( naold != na ) {
        naold = na;
        kloold = 1;
        khiold = na;
    }

    klo = ( x > xa[kloold] ? kloold : 1 );
    khi = ( x < xa[khiold] ? khiold : na );
    while ( khi - klo > 1) {

        k = (khi + klo) >> 1;
        if ( xa[k] > x ) khi = k;
        else    klo = k;
    }
    kloold = klo;
    khiold = khi;

    h = xa[khi] - xa[klo];
    if ( h == 0.0 ) nrerror( "Bad xa input to routine spline_interp().");
    a = (xa[khi] - x)/h;
    b = (x - xa[klo])/h;
    return a*ya[klo] + b*ya[khi] + (a*(a*a - 1.0)*yapp[klo] + b*(b*b - 1.0)*yapp[khi])*(h*h)/6.0;
}




/**************************************************************************************************************************

                linear_interp

        Given the arrays xa[1..na], ya[1..na] which tabulate a function (with the xa's in order), 
        and given a value of x, this routine returns a simple linear interpolated value y.

**************************************************************************************************************************/


double linear_interp( double *xa, double *ya, int na, double x ) {

    int k, klo, khi;
    double h, b, a;

    /*
     *    Store the last values of klo and khi to see if they still work.
     *    Also store the last number of points to watch out for new spline data.
     *    If the old points don't work reset them to the endpoints.
     *    In either case, refine klo and khi with bisection.
     */
    static int naold = 0;
    static int kloold;
    static int khiold;

    /*
     *    If the number of points has changed then spline was rerun or a new set of points
     *    is being provided.  To prevent the arrays being indexed to an illegal value and
     *    to reflect the arrival of new data reset the values of kloold and khiold and naold.
     */
    if ( naold != na ) {
        naold = na;
        kloold = 1;
        khiold = na;
    }

    klo = ( x > xa[kloold] ? kloold : 1 );
    khi = ( x < xa[khiold] ? khiold : na );
    while ( khi - klo > 1) {

        k = (khi + klo) >> 1;
        if ( xa[k] > x) khi = k;
        else    klo = k;
    }
    kloold = klo;
    khiold = khi;

    h = xa[khi] - xa[klo];
    if ( h == 0.0 ) nrerror( "Bad xa input to routine linear_interp().");
    a = (xa[khi] - x)/h;
    b = (x - xa[klo])/h;
    return a*ya[klo] + b*ya[khi];
}




