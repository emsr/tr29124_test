
#include <math.h>


#include "nric.h"




/*******************************************************************************

    This routine returns the spherical Bessel function of order 0.

*******************************************************************************/

double sph_bessel_j0( double x ) {

    double  j0;

    if  ( x < 0.0 ) nrerror("Negative argument in sph_bessel_j0." );

    if  ( x == 0.0 ) return 1.0;

    j0 = sin(x)/x;
    return  j0;
}



/*****************************************************************************************************

    This routine returns the spherical Bessel function of order 1.

*****************************************************************************************************/

double sph_bessel_j1( double x ) {

    double  j1;

    if  ( x < 0.0 ) nrerror( "Negative argument in sph_bessel_j1." );

    if  ( x == 0.0 ) return 0.0;

    j1 = ( sin(x) - x*cos(x) )/(x*x);
    return  j1;
}




/*****************************************************************************************************

    This routine returns the spherical Bessel function of order n.

*****************************************************************************************************/

double  sph_bessel_j( int n, double x) {

    int  m;
    double  j0, j1, ans, jp1, j, jm1, sum;

    const double SMALL = 1.0e-6;
    const double LARGE = 1.0e+6;

    const double ACC = 256.0;
    const double BIGNO = 1.0e10;
    const double BIGNI = 1.0e-10;


    if  ( n < 0 ) nrerror( "Negative order in sph_bessel_j." );
    if  ( x < 0.0 ) nrerror( "Negative argument in sph_bessel_j." );
    if  ( ( x == 0.0 ) && ( n == 0 ) ) return 1.0;
    if  ( ( x == 0.0 ) && ( n != 0 ) ) return 0.0;


    /*
     *    small argument value
     */
    if  ( x < SMALL )    {

        ans = 1.0;
        for ( m = 1; m <= n; m++ ) ans *= x/(2.0*m + 1.0);
        return ans;
    }

    /*
     *    large argument value
     */
    if  ( x > LARGE*n*(n + 1) )  return sin( x - 0.5*PI*n )/x;

    /*
     *    Compute j0
     */
    j0 = sin(x)/x;
    if  ( n == 0 ) return j0;

    /*
     *    Compute j1
     */
    j1 = ( sin(x) - x*cos(x) )/(x*x);
    if  ( n == 1) return j1;

    /*
     *    Compute ans, n >= 2
     */
    if ( fabs(x) > (1.0*n + 0.5) ) {

        /*
         *    Use forward recursion.
         */
        for  ( ans = 0.0, jm1 = j0, j = j1, m = 1; m < n; ++m )    {

            jp1 = ( 2*m + 1 )*j/x - jm1;
            jm1 = j;
            j = jp1;
        }
        ans = j;

    } else {

        /*
         *    Use Miller's backward recursion.
         */
        for ( sum = 0.0, ans = 0.0, j = 1.0, jp1 = 0.0, m = n + (int) sqrt(ACC*n); m >= 1; --m ) {

            jm1 = ( 2*m + 1 )*j/x - jp1;
            jp1 = j;
            j = jm1;
            if ( fabs(j) > BIGNO ) {
                /*
                 *    Renormalize to prevent overflows.
                 */
                j *= BIGNI;
                jp1 *= BIGNI;
                ans *= BIGNI;
                sum *= BIGNI*BIGNI;
            }
            sum += (2*m - 1)*j*j;
            if ( m == n )  ans = jp1;
        }
        ans /= sqrt(sum);
    }
    return  ans;
}




/*****************************************************************************************************

    This routine returns the spherical Bessel function of order 0.

*****************************************************************************************************/

double sph_bessel_y0( double x ) {

    double  y0;

    if  ( x < 0.0 )  nrerror( "Negative argument in sph_bessel_y0." );

    if  ( x == 0.0 )  return -1.0;

    y0 = -cos(x)/x;
    return  y0;
}




/*****************************************************************************************************

    This routine returns the spherical Bessel function of order 1.

*****************************************************************************************************/

double sph_bessel_y1( double x ) {

    double y1;

    if  ( x < 0.0 ) nrerror( "Negative argument in sph_bessel_y1." );

    if  ( x == 0.0 ) return 1.0;

    y1 = -( cos(x) + x*sin(x) )/(x*x);
    return  y1;
}




/*****************************************************************************************************

    This routine returns the spherical Bessel function of order n.

*****************************************************************************************************/

double sph_bessel_y( int n, double x) {

    int  m;
    double  y0, y1, sum, ans, y, ym1, yp1;
    const double SMALL = 1.0e-6;
    const double LARGE = 1.0e+6;

    if  ( n < 0 ) nrerror( "Negative order in sph_bessel_y." );
    if  ( x < 0.0 ) nrerror( "Negative argument in sph_bessel_y." );
    if  ( ( x == 0.0 ) && ( n == 0 ) ) return 1.0;
    if  ( ( x == 0.0 ) && ( n != 0 ) ) return 0.0;


    /*
     *    small argument value
     */
    if  ( x < SMALL )    {

        ans = -1.0/x;
        for ( m = 1; m <= n; m++ )  ans *= (2*m - 1)/x;
        return ans;
    }

    /*
     *    large argument value
     */
    if  ( x > LARGE*n*(n + 1) )  return  -cos(x - 0.5*PI*n)/x;

    /*
     *    Compute y0
     */
    y0 = -cos(x)/x;
    if  ( n == 0 )  return y0;

    /*
     *    Compute y1
     */
    y1 = -( cos(x) + x*sin(x) )/(x*x);
    if  ( n == 1 )  return y1;

    /*
     *    Compute ans, n >= 2
     */
    for  ( ans = 0.0, ym1 = y0, y = y1, m = 1; m < n; m++ )    {

        yp1 = (2*m + 1)*y/x - ym1;
        ym1 = y;
        y = yp1;
    }
    ans = y;
    return  ans;
}




/*******************************************************************************

    This routine returns the associated spherical Bessel function of order 0.

*******************************************************************************/

double sph_bessel_i0( double x ) {

    double  i0;

    if  ( x < 0.0 ) nrerror("Negative argument in sph_bessel_i0." );

    if  ( x == 0.0 ) return 1.0;

    i0 = sinh(x)/x;
    return  i0;
}



/*****************************************************************************************************

    This routine returns the associated spherical Bessel function of order 1.

*****************************************************************************************************/

double sph_bessel_i1( double x ) {

    double  i1;

    if  ( x < 0.0 ) nrerror( "Negative argument in sph_bessel_i1." );

    if  ( x == 0.0 ) return 0.0;

    i1 = ( -sinh(x) + x*cosh(x) )/(x*x);
    return  i1;
}



/*****************************************************************************************************

    This routine returns the associated spherical Bessel function of order n.

*****************************************************************************************************/

double sph_bessel_i( int n, double x ) {

    int  m;
    double  i0, i1, sum, ans, i, im1, ip1;

    const double ACC = 256.0;
    const double BIGNO = 1.0e10;
    const double BIGNI = 1.0e-10;

    if  ( n < 0 ) nrerror( "Negative order in sph_bessel_i." );
    if  ( x < 0.0 ) nrerror( "Negative argument in sph_bessel_i." );
    if  ( ( x == 0.0 ) && ( n == 0 ) ) return 1.0;
    if  ( ( x == 0.0 ) && ( n != 0 ) ) return 0.0;

    /*
     *    Compute i0
     */
    i0 = sinh(x)/x;
    if  ( n == 0 )  return i0;

    /*
     *    Compute i1
     */
    i1 = ( -sinh(x) + x*cosh(x) )/(x*x);
    if  ( n == 1 )  return i1;

    /*
     *    Compute ans, n >= 2
     *    Use Miller's backward recursion.
     */
    if ( fabs(x) > 1.0*n + 0.5 ) {
        for  ( ans = 0.0, im1 = i0, i = i1, m = 1; m < n; m++ )    {

            ip1 = (2*m + 1)*i/x + im1;
            im1 = i;
            i = ip1;
        }
        ans = i;
    } else {
        for ( sum = 0.0, ans = 0.0, i = 1.0, ip1 = 0.0, m = n + (int) sqrt(ACC*n); m >= 1; --m ) {

            im1 = ( 2*m + 1 )*i/x + ip1;
            ip1 = i;
            i = im1;
            if ( fabs(i) > BIGNO ) {
                /*
                 *    Renormalize to prevent overflows.
                 */
                i *= BIGNI;
                ip1 *= BIGNI;
                ans *= BIGNI;
                sum *= BIGNI*BIGNI;
            }
            sum += (2*m - 1)*i*i;
            if ( m == n )  ans = ip1;
        }
        ans /= sqrt(sum);
    }

    return ans;
}




/*******************************************************************************

    This routine returns the associated spherical Bessel function of order 0.

*******************************************************************************/

double sph_bessel_k0( double x ) {

    double  k0;

    if  ( x < 0.0 ) nrerror("Negative argument in sph_bessel_k0." );

    if  ( x == 0.0 ) return 1.0;

    k0 = (PI/2.0)*exp(-x)/x;
    return  k0;
}



/*****************************************************************************************************

    This routine returns the associated spherical Bessel function of order 1.

*****************************************************************************************************/

double sph_bessel_k1( double x ) {

    double  k1;

    if  ( x < 0.0 ) nrerror( "Negative argument in sph_bessel_k1." );

    if  ( x == 0.0 ) return 0.0;

    k1 = (PI/2.0)*exp(-x)*( 1.0 + 1.0/x )/x;
    return  k1;
}



/*****************************************************************************************************

    This routine returns the associated spherical Bessel function of order n.

*****************************************************************************************************/

double sph_bessel_k( int n, double x ) {

    int m;
    double k0, k1, km1, k, kp1, sum, ans;

    if  ( n < 0 ) nrerror( "Negative order in sph_bessel_k." );
    if  ( x < 0.0 ) nrerror( "Negative argument in sph_bessel_k." );
    if  ( ( x == 0.0 ) && ( n == 0 ) ) return 1.0;
    if  ( ( x == 0.0 ) && ( n != 0 ) ) return 0.0;

    /*
     *    Compute k0
     */
    k0 = (PI/2.0)*exp(-x)/x;
    if  ( n == 0 )  return k0;

    /*
     *    Compute k1
     */
    k1 = (PI/2.0)*exp(-x)*( 1.0 + 1.0/x )/x;
    if  ( n == 1 )  return k1;

    /*
     *    Compute ans, n >= 2
     */
    for  ( ans = 0.0, km1 = k0, k = k1, m = 1; m < n; m++ )    {

        kp1 = (2*m + 1)*k/x + km1;
        km1 = k;
        k = kp1;
    }
    ans = k;

    return k;
}


