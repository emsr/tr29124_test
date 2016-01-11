


#include <math.h>


#include "nric.h"


/*
 *    Chebyshev fit.  Given a function, the lower limit a, the upper limit b,
 *    the number of points n this routine computes the coefficients  c[0..n-1] of a 
 *    Chebyshev polynomial expansion such that funk(x) =~ kSum_0^n-1 c_k T_k(x) - c_0/2.
 *    Use this routine with moderately large n of 30 or 50.  Then you can truncate the
 *    series to a smaller number of terms to satisfy the accuracy requirements.
 */
void chebyshev_fit( double a, double b, double *c, int n, double (*funk)( double ) ) {

    int j, k;
    double fac, bma, bpa, *f;

    f = dvector( 0, n-1 );

    bma = 0.5*(b - a);
    bpa = 0.5*(b + a);
    for ( k = 0; k < n; ++k ) {
        double y = cos( PI*( k + 0.5 )/n );
        f[k] = (*funk)( bpa + y*bma );
    }

    fac = 2.0/n;
    for ( j = 0; j < n; ++j ) {
        double sum = 0.0;
        for ( k = 0; k < n; ++k ) sum += f[k]*cos( PI*j*( k + 0.5 )/n );
        c[j] = fac*sum;
    }

    free_dvector( f, 0, n-1 );
}


/*
 *    Chebyshev evaluation.  All arguments are input.  c[0..m-1] is an array of Chebyshev
 *    coefficients - the first m elements of c output from chebyshev_fit (previously called
 *    with the same a and b).  The Chebyshev polynomial is evaluated at a point y determined
 *    from x, a, b.  The result is returned as the function value.
 */
double chebyshev_eval( double a, double b, double *c, int m, double x ) {

    int j;
    double d = 0.0, dd = 0.0, sv, y, y2;

    if ( ( x - a )*( x - b ) > 0.0 ) nrerror( "x not in range in routine chebyshev_eval." );

    y2 = 2.0*( y = ( 2.0*x - a - b )/( b - a ) );
    for ( j = m - 1; j >= 1; --j ) {
        sv = d;
        d = y2*d - dd + c[j];
        dd = sv;
    }
    return y*d - dd + 0.5*c[0];
}



/*
 *    Given a, b, and c as output from chebyshev_fit, and given m, the desired degree of
 *    approximation, this routine returns the array cder[0..m-1], the Chebyshev coefficients
 *    of the derivative of the function whose coefficients are c.
 */
void chebyshev_deriv( double a, double b, double *c, double *cder, int m ) {

    int j;
    double con;

    cder[m-1] = 0.0;
    cder[m-2] = 2*(m - 1)*c[m-1];
    for ( j = m-3; j >= 0; --j ) cder[j] = cder[j+2] + 2*(j + 1)*c[j+1];
    con = 2.0/(b - a);
    for ( j = 0; j < m; ++j ) cder[j] *= con;
}



/*
 *    Given a, b, and c as output from chebyshev_fit, and given m, the desired degree of
 *    approximation, this routine returns the array cint[0..m-1], the Chebyshev coefficients
 *    of the integral of the function whose coefficients are c.
 *    The constant of integration is set so that the integral vanishes at a.
 */
void chebyshev_integ( double a, double b, double *c, double *cint, int m ) {

    int j;
    double sum = 0.0, fac = 1.0, con;

    con = 0.25*(b - a);
    for ( j = 1; j < m-1; ++j ) {
        /*
         *    Accumulate the constant of integration in sum.
         *    fac will equal +/-1 to alternate the sign of the terms.
         */
        cint[j] = con*(c[j-1] - c[j+1])/j;
        sum += fac*cint[j];
        fac = -fac;
    }
    cint[m-1] = con*c[m-2]/(m - 1);
    sum += fac*cint[m-1];
    /*
     *    Set the constant of integration.
     */
    cint[0] = 2.0*sum;
}



/*
 *    Given a, b, and c as output from chebyshev_fit, and given m, the desired degree of
 *    approximation, this routine returns the array d[0..m-1], of coefficients of a polynomial
 *    expansion which is equivalent to the Chebyshev fit.
 */
void chebyshev_2_poly( double a, double b, double *c, double *d, int m ) {

    int j, k;
    double sv, *dd;

    dd = dvector( 0, m-1 );

    for ( j = 0; j < m; ++j ) d[j] = dd[j] = 0.0;
    d[0] = c[m-1];
    for ( j = m-2; j >= 1; --j ) {
        for ( k = m-j; k >= 1; --k ) {
            sv = d[k];
            d[k] = 2.0*d[k-1] - dd[k];
            dd[k] = sv;
        }
        sv = d[0];
        d[0] = -dd[0] + c[j];
        dd[0] = sv;
    }
    for ( j = m-1; j >= 1; --j ) d[j] = d[j-1] - dd[j];
    d[0] = -dd[0] + 0.5*c[0];

    free_dvector( dd, 0, m-1 );

    /*
     *    Map the interval [-1,+1] to [a,b].
     */
    poly_shift_coeff( a, b, d, m );
}


/*
 *    Converts a polynomial to a Chebyshev expansion.
 *    A polynomial valid in the range x=a to x=b is specified by coefficients d[0..m-1].
 *    An equivalent array of Chebyshev coefficients, c[0..m-1], is output.
 *    The index mfew will be set to the index of the first nonzero Chebyshev coefficient smaller than err.
 *    Then the polynomial is approximated by the first mfew cefficients c[0..mfew-1].
 */
void poly_2_chebyshev( double a, double b, double *d, double *c, int m ) {

    int j, jm, jp, k;
    double fac, pow;

    if ( b == a ) nrerror( "Input range must be nonzero in poly_2_chebyshev." );

    /*
     *    Map the interval of the polynomial from x:[a,b] to y:[-1,+1].
     */
    poly_shift_coeff( (-2.0 - b - a)/(b - a), (+2.0 - b - a)/(b - a), d, m );

    pow = 1.0;
    c[0] = 2.0*d[0];
    for ( k = 1; k < m; ++k ) {
        c[k] = 0.0;
        fac = d[k]/pow;
        jm = k;
        jp = 1;
        for ( j = k; j >= 0; j-=2, --jm, ++jp ) {
            /*
             *    Increment this and lower orders of Chebyshev with the combinatorial
             *    coefficient times d[k].
             */
            c[j] += fac;
            fac *= ((double)jm)/((double)jp);
        }
        pow += pow;
    }

    /*
     *    Map the interval of the polynomial from y:[-1,+1] to x:[a,b].
     */
    poly_shift_coeff( a, b, d, m );
}


/*
 *    Economize the polynomial d[0..m-1] by a new polynomial d[0..m-1] so that the new polynomial
 *    evaluated with fewer terms d[0..mfew-1] will equal the old polynomial within the specified tolerance err.
 *    The index mfew will be set to the index of the first nonzero Chebyshev coefficient smaller than err.
 *    Then the polynomial is approximated by the first mfew cefficients d[0..mfew-1].
 */
void poly_economize( double a, double b, double *d, double *c, int m, int *mfew, double err ) {

    int k;

    /*
     *    Convert the shifted polynomial into a Chebyshev series.
     */
    poly_2_chebyshev( a, b, d, c, m );

    /*
     *    Truncate the Chebyshev series so that the first nonzero neglected term is less than err.
     *    Skipping zero coefficients takes care of series with alternating 0 and nonzero terms.
     *    Truncate the polynomial as well.
     */
    for ( *mfew = 0, k = 0; k < m; ++k ) {
        if ( c[k] == 0.0 ) continue;
        else if ( c[k] < err ) {
            *mfew = k;
            break;
        } else *mfew = k+1;
    }
    for ( k = *mfew; k < m; ++k ) d[k] = c[k] = 0.0;

    /*
     *    Convert the truncated Chebyshev series back into a polynomial.
     */
    chebyshev_2_poly( a, b, c, d, *mfew );
}


/*
 *    Evaluate the definite integral of a function fitted by Chebyshev polynomials
 *    over the interval from a to b.  The output coefficients from chebyshev_fit
 *    are in c.  The maximum order used for the integral is m.
 */
double clenshaw_curtis_quad( double a, double b, double *c, int m, double eps ) {

    int k;
    double sum;

    if ( eps <= 0.0 ) nrerror( "Error tolerance eps must be greater than 0 in clenshaw_curtis_quad." );

    for ( sum = 0.0, k = 0; (2*k + 1) < m; ++k ) sum -= c[2*k + 1]/((2*k + 1)*(2*k - 1));
    sum -= c[1]/2.0;

    return (b - a)*sum;
}


