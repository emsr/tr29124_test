

#include <math.h>

#include "nric.h"


/*
 *    Returns the Bessel function J of order 0.
 *    The accuracy is really 10e-7 or so or single precision.
 */
double bessel_j0( double x )  {

    double ax, z;
    double xx, y, ans, ans1, ans2;

    if  ( ( ax = fabs( x ) ) < 8.0 ) {
        y = x*x;
        ans1 = 57568490574.0 + y*( -13362590354.0 + y*( 651619640.7 + y*( -11214424.18 + y*( 77392.33017 + y*( -184.9052456 ) ) ) ) );
        ans2 = 57568490411.0 + y*( 1029532985.0 + y*( 9494680.718 + y*( 59272.64853 + y*( 267.8532712 + y ) ) ) );
        ans = ans1/ans2;
    } else {
        z = 8.0/ax;
        y = z*z;
        xx = ax - 0.785398164;
        ans1 = 1.0 + y*( -0.1098628627e-2 + y*( 0.2734510407e-4 + y*( -0.2073370639e-5 + y*( 0.2093887211e-6 ) ) ) );
        ans2 = -0.1562499995e-1 + y*( 0.1430488765e-3 + y*( -0.6911147651e-5 + y*( 0.7621095161e-6 + y*( -0.934935152e-7 ) ) ) );
        ans = sqrt(0.636619772/ax)*( cos(xx)*ans1 - z*sin(xx)*ans2 );
    }
    return ans;
}



/*
 *    Returns the Bessel function J of order 1.
 *    The accuracy is really 10e-7 or so or single precision.
 */
double bessel_j1( double x )  {

    double ax, z;
    double xx, y, ans, ans1, ans2;

    if  ( ( ax = fabs( x ) ) < 8.0 ) {
        y = x*x;
        ans1 = x*( 72362614232.0 + y*( -7895059235.0 + y*( 242396853.1
              + y*( -2972611.439 + y*( 15704.48260 + y*( -30.16036606 ) ) ) ) ) );
        ans2 = 144725228442.0 + y*( 2300535178.0 + y*( 18583304.74 + y*( 99447.43394 + y*( 376.9991397 + y ) ) ) );
        ans = ans1/ans2;
    } else {
        z = 8.0/ax;
        y = z*z;
        xx = ax - 2.356194491;
        ans1 = 1.0 + y*( 0.183105e-2 + y*( -0.3516396496e-4 + y*( 0.2457520174e-5 + y*( -0.240337019e-6 ) ) ) );
        ans2 = 0.4687499995e-1 + y*( -0.2002690873e-3 + y*( 0.8449199096e-5 + y*( -0.88228987e-6 + y*( 0.105787412e-6 ) ) ) );
        ans = sqrt(0.636619772/ax)*( cos(xx)*ans1 - z*sin(xx)*ans2 );
    }
    return ans;
}



/*
 *    Returns the Bessel function J of order n.
 *    The accuracy is really 10e-7 or so or single precision.
 */
double bessel_j( int n, double x ) {

    int j, jsum, m;
    double ax, bj, bjm, bjp, sum, tox, ans;

    const double ACC = 256.0;
    const double BIGNO = 1.0e10;
    const double BIGNI = 1.0e-10;

    if ( n < 0 ) nrerror( "Negative order in bessel_j." );

    if ( n == 0 ) return bessel_j0(x);

    if ( n == 1 ) return bessel_j1(x);

    ax = fabs(x);
    if ( ax == 0.0 ) {
        return 0.0;
    } else if ( ax > (double) n ) {
        /*
         *    Use the forward recursion if |n/x| < 1.
         */
        tox = 2.0/x;
        bjm = bessel_j0(ax);
        bj = bessel_j1(ax);
        for ( j = 1; j < n; ++j ) {
            bjp = j*tox*bj - bjm;
            bjm = bj;
            bj = bjp;
        }
        ans = bj;
    } else {
        /*
         *    Otherwise, use Miller's backward recursion.
         */
        tox = 2.0/x;
        /*
         *    This multiplying and dividing by 2 makes m an even number.
         */
        m = 2*( n + (int) sqrt(ACC*n) )/2;
        jsum = 0;
        bjp = ans = sum = 0.0;
        bj = 1.0;
        for ( j = m; j > 0; --j ) {
            bjm = j*tox*bj - bjp;
            bjp = bj;
            bj = bjm;
            if ( fabs(bj) > BIGNO ) {
                /*
                 *    Renormalize to prevent overflows.
                 */
                bj *= BIGNI;
                bjp *= BIGNI;
                ans *= BIGNI;
                sum *= BIGNI;
            }
            if ( jsum ) sum += bj;
            jsum = !jsum;
            if ( j == n ) ans = bjp;
        }
        sum = 2.0*sum - bj;
        ans /= sum;
    }
    return (x < 0.0) && n%2 == 1 ? -ans : ans;
}



/*
 *    Returns the Bessel function Y of order 0.
 *    The accuracy is really 10e-7 or so or single precision.
 */
double bessel_y0( double x )  {

    double z;
    double xx, y, ans, ans1, ans2;

    if ( x <= 0.0 ) nrerror( "Negative argument in bessel_y0." );

    if  ( x < 8.0 ) {
        y = x*x;
        ans1 = -2957821389.0 + y*( 7062834065.0 + y*( -512359803.6 + y*( 10879881.29 + y*( -86327.92757 + y*( 228.4622733 ) ) ) ) );
        ans2 = 40076544269.0 + y*( 745249964.8 + y*( 7189466.438 + y*( 47447.26470 + y*( 226.1030244 + y ) ) ) );
        ans = ans1/ans2 + 0.636619772*bessel_j0(x)*log(x);
    } else {
        z = 8.0/x;
        y = z*z;
        xx = x - 0.785398164;
        ans1 = 1.0 + y*( -0.1098628627e-2 + y*( 0.2734510407e-4 + y*( -0.2073370639e-5 + y*( 0.2093887211e-6 ) ) ) );
        ans2 = -0.1562499995e-1 + y*( 0.1430488765e-3 + y*( -0.6911147651e-5 + y*( 0.7621095161e-6 + y*( -0.934945152e-7 ) ) ) );
        ans = sqrt(0.636619772/x)*( sin(xx)*ans1 + z*cos(xx)*ans2 );
    }
    return ans;
}




/*
 *    Returns the Bessel function Y of order 1.
 *    The accuracy is really 10e-7 or so or single precision.
 */
double bessel_y1( double x )  {

    double z;
    double xx, y, ans, ans1, ans2;

    if ( x <= 0.0 ) nrerror( "Negative argument in bessel_y1." );

    if  ( x < 8.0 ) {
        y = x*x;
        ans1 = x*( -0.4900604943e13 + y*( 0.1275274390e13 + y*( -0.5153438139e11 + y*( 0.7349264551e9
               + y*( -0.4237922726e7 + y*( 0.8511937935e4 ) ) ) ) ) );
        ans2 = 0.2499580570e14 + y*( 0.4244419664e12 + y*( 0.3733650367e10 + y*( 0.2245904002e8
               + y*( 0.1020426050e6 + y*( 0.3549632885e3 + y ) ) ) ) );
        ans = ans1/ans2 + 0.636619772*( bessel_j1(x)*log(x) - 1.0/x );
    } else {
        z = 8.0/x;
        y = z*z;
        xx = x - 2.356194491;
        ans1 = 1.0 + y*( 0.183105e-2 + y*( -0.3516396496e-4 + y*( 0.2457520174e-5 + y*( -0.240337019e-6 ) ) ) );
        ans2 = 0.4687499995e-1 + y*( -0.2002690873e-3 + y*( 0.8449199096e-5 + y*( -0.88228987e-6 + y*( 0.105787412e-6 ) ) ) );
        ans = sqrt(0.636619772/x)*( sin(xx)*ans1 + z*cos(xx)*ans2 );
    }
    return ans;
}



/*
 *    Returns the Bessel function Y of order n.
 *    The accuracy is really 10e-7 or so or single precision.
 */
double bessel_y( int n, double x ) {

    int j;
    double y0, y1, yp, tox;

    if ( x <= 0.0 ) nrerror( "Negative argument in bessel_y1." );

    if ( n < 0 ) nrerror( "Negative order in bessel_y." );

    y0 = bessel_y0(x);
    if ( n == 0 ) return y0;

    y1 = bessel_y1(x);
    if ( n == 1 ) return y1;

    tox = 2.0/x;
    for( j = 1; j <= n; ++j ) {
        yp = j*tox*y1 - y0;
        y0 = y1;
        y1 = yp;
    }
    return y1;
}



/*
 *    Returns the associated Bessel function I of order 0.
 *    The accuracy is really 10e-7 or so or single precision.
 */
double bessel_i0( double x ) {

    double ax, ans;
    double y;

    if ( ( ax = fabs(x) ) < 3.75 ) {
        y = x/3.75;
        y *= y;
        ans = 1.0 + y*( 3.5156229 + y*( 3.0899424 + y*( 1.2067492 + y*( 0.2659732 + y*( 0.360768e-1 + y*( 0.45813e-2 ) ) ) ) ) );
    } else {
        y = 3.75/ax;
        ans = 0.39894228 + y*( 0.1328592e-1 + y*( 0.225319e-2 + y*( -0.157565e-2 + y*( 0.916281e-2
                         + y*( -0.2057706e-1 + y*( 0.2635537e-1 + y*( -0.1647633e-1 + y*( 0.392337e-2 ) ) ) ) ) ) ) );
        ans *= exp(ax)/sqrt(ax);
    }
    return ans;
}



/*
 *    Returns the associated Bessel function I of order 1.
 *    The accuracy is really 10e-7 or so or single precision.
 */
double bessel_i1( double x ) {

    double ax, ans;
    double y;

    if ( ( ax = fabs(x) ) < 3.75 ) {
        y = x/3.75;
        y *= y;
        ans = ax*( 0.5 + y*( 0.87890594 + y*( 0.51498869 + y*( 0.15084934
             + y*( 0.2658733e-1 + y*( 0.301532e-2 + y*( 0.32411e-3 ) ) ) ) ) ) );
    } else {
        y = 3.75/ax;
        ans = 0.39894228 + y*( -0.3988024e-1 + y*( -0.362018e-2 + y*( 0.163801e-2 + y*( -0.1031555e-1
                         + y*( 0.2282967e-1 + y*( -0.2895312e-1 + y*( 0.1787654e-1 + y*( -0.420059e-2 ) ) ) ) ) ) ) );
        ans *= exp(ax)/sqrt(ax);
    }
    return x < 0.0 ? -ans : ans;
}




/*
 *    Returns the associated Bessel function I of order n.
 *    The accuracy is really 10e-7 or so or single precision.
 */
double bessel_i( int n, double x ) {

    int j, m;
    double bi, bim, bip, tox, ans;

    const double ACC = 256.0;
    const double BIGNO = 1.0e10;
    const double BIGNI = 1.0e-10;

    if ( n < 0 ) nrerror( "Negative order in bessel_i." );

    if ( n == 0 ) return bessel_i0(x);

    if ( n == 1 ) return bessel_i1(x);

    if ( x == 0.0 ) {
        return 0.0;
    } else if ( fabs(x) > (double) n ) {
        /*
         *    Use forward recurrence.
         */
        tox = 2.0/x;
        for ( j = 1; j < n; ++j ) {
            bip = j*tox*bi + bim;
            bim = bi;
            bi = bip;
        }
        ans = bip;
    } else {
        /*
         *    Use Miller's backward recursion to compute I_m through I_0 (m > n).
         *    Save I_n but we need I_0 to renormalize the final answer I_n.
         */
        m = 2*( n + (int) sqrt(ACC*n) )/2;
        tox = 2.0/fabs(x);
        bip = ans = 0.0;
        bi = 1.0;
        for ( j = m; j > 0; --j ) {
            bim = j*tox*bi + bip;
            bip = bi;
            bi = bim;
            if ( fabs(bi) > BIGNO ) {
                /*
                 *    Renormalize to prevent overflows.
                 */
                bi *= BIGNI;
                bip *= BIGNI;
                ans *= BIGNI;
            }
            if ( j == n ) ans = bip;
        }
        ans *= bessel_i0(x)/bi;
    }
    return x < 0.0 && n%2 == 1 ? -ans : ans;
}



/*
 *    Returns the associated Bessel function I of order 0.
 *    The accuracy is really 10e-7 or so or single precision.
 */
double bessel_k0( double x ) {

    double y, ans;

    if ( x <= 0.0 ) nrerror( "Negative argument in bessel_k0." );

    if ( x <= 2.0 ) {
        y = x*x/4.0;
        ans = -log(x/2.0)*bessel_i0(x) + -0.57721566 + y*( 0.42278420 + y*( 0.23069756 + y*( 0.3488590e-1
              + y*( 0.262698e-2 + y*( 0.10750e-3 + y*( 0.74e-5 ) ) ) ) ) ) ;
    } else {
        y = 2.0/x;
        ans = 1.25331414 + y*( -0.7832358e-1 + y*( 0.2189568e-1 + y*( -0.1062446e-1 
              + y*( 0.587872e-2 + y*( -0.251540e-2 + y*( 0.53208e-3 ) ) ) ) ) );
        ans *= exp(-x)/sqrt(x);
    }
    return ans;
}



/*
 *    Returns the associated Bessel function I of order 1.
 *    The accuracy is really 10e-7 or so or single precision.
 */
double bessel_k1( double x ) {

    double y, ans;

    if ( x <= 0.0 ) nrerror( "Negative argument in bessel_k1." );

    if ( x <= 2.0 ) {
        y = x*x/4.0;
        ans = log(x/2.0)*bessel_i1(x) + (1.0/x)*( 1.0 + y*( 0.15443144 + y*( -0.67278579 + y*( -0.18156897
              + y*( -0.1919402e-1 + y*( -0.110404e-2 + y*( -0.4686e-4 ) ) ) ) ) ) ) ;
    } else {
        y = 2.0/x;
        ans = 1.25331414 + y*( 0.23498619 + y*( -0.3655620e-1 + y*( 0.1504268e-1 
              + y*( -0.780353e-2 + y*( 0.325614e-2 + y*( -0.68245e-3 ) ) ) ) ) );
        ans *= exp(-x)/sqrt(x);
    }
    return ans;
}



/*
 *    Returns the associated Bessel function I of order n.
 *    The accuracy is really 10e-7 or so or single precision.
 */
double bessel_k( int n, double x ) {

    int j;
    double bk, bkm, bkp, tox;

    if ( x <= 0.0 ) nrerror( "Negative argument in bessel_k." );

    if ( n < 0 ) nrerror( "Negative order in bessel_k." );

    bkm = bessel_k0(x);
    if ( n == 0 ) return bkm;

    bk = bessel_k1(x);
    if ( n == 1 ) return bk;

    tox = 2.0/x;
    for ( j = 1; j < n; ++j ) {
        bkp = j*tox*bk + bkm;
        bkm = bk;
        bk = bkp;
    }
    return bk;
}


