

#include <math.h>

#include "nric.h"


static double Sin( double x ) { return sin( PIO2*x*x ); }

static double Cos( double x ) { return cos( PIO2*x*x ); }


/*
 *    Returns the Fresnel cosine integral by straight Romberg integration.
 */
double fresnel_c( double x ) { return quad_romberg( Cos, 0.000000001, x, 0.0000001 ); }


/*
 *    Returns the Fresnel sine integral by straight Romberg integration.
 */
double fresnel_s( double x ) { return quad_romberg( Sin, 0.000000001, x, 0.0000001 ); }


/*
 *    Computes the Fresnel integrals C(x) and S(x) for all x by either series expansion or
 *    by complex continued fraction as appropriate.
 */
double fresnel( double x, double *c, double *s ) {

    int k, n, odd;
    double a, ax, fact, pix2, sign, sum, sums, sumc, term, test;
    dcomplex b, cc, d, h, del, cs;

    const double EPS = 6.0e-8;
    const int MAXIT = 100;
    const double FPMIN = 1.0e-30;
    const double XMIN = 1.5;
    const int TRUE = 1;

    const dcomplex ONE = mk_dc( 1.0, 0.0 );
    const dcomplex FOUR = mk_dc( 4.0, 0.0 );

    ax = fabs(x);
    if ( ax < sqrt(FPMIN) ) {
        *s = 0.0;
        *c = ax;
    }  else if ( ax < XMIN ) {
        /*
         *    Evaluate S and C by series expansion.
         */
        sum = sums = 0.0;
        sumc = ax;
        sign = 1.0;
        fact = PIO2*ax*ax;
        odd = TRUE;
        term = ax;
        n = 3;
        for ( k = 1; k <= MAXIT; ++k ) {
            term *= fact/k;
            sum += sign*term/n;
            test = fabs(sum)*EPS;
            if ( odd ) {
                sign = -sign;
                sums = sum;
                sum = sumc;
            } else {
                sumc = sum;
                sum = sums;
            }
            if ( term < test ) break;
            odd = !odd;
            n += 2;
        }
        if ( k > MAXIT ) nrerror( "Series evaluation failed in fresnel." );
        *c = sumc;
        *s = sums;
    } else {
        /*
         *    Evaluate S and C by Lentz's complex continued fraction method.
         */
        pix2 = PI*ax*ax;
        b = mk_dc( 1.0, -pix2 );
        cc = mk_dc( 1.0/FPMIN, 0.0 );
        d = h = div_dc( ONE, b );
        n = -1;
        for ( k = 2; k <= MAXIT; ++k ) {
            n += 2;
            a = -n*(n+1);
            b = add_dc( b, FOUR );
            d = div_dc( ONE, add_dc( rmul_dc( a, d ), b ) );
            cc = add_dc( b, div_dc( mk_dc( a, 0.0 ), cc ) );
            del = mul_dc( cc, d );
            h = mul_dc( h, del );
            if ( fabs( del.r - 1.0 ) + fabs( del.i ) < EPS ) break;
        }
        if ( k > MAXIT ) nrerror( "Continued fraction evaluation failed in fresnel." );
        h = mul_dc( mk_dc( ax, -ax ), h );
        cs = mul_dc( mk_dc( 0.5, 0.5 ), sub_dc( ONE, mul_dc( mk_dc( cos(0.5*pix2), sin(0.5*pix2) ), h ) ) );
        *c = cs.r;
        *s = cs.i;
    }
    if ( x < 0.0 ) {
        *c = -*c;
        *s = -*s;
    }

    return 0.0;
}



double fock_r( double x ) {

    double sqrtx, c, s;

    sqrtx = sqrt(x);

    fresnel( sqrtx/SQRTPIO2, &c, &s );

    /*
    s = quad_romberg( Sin, 0.000000001, sqrtx/SQRTPIO2, 0.0000001 );
    c = quad_romberg( Cos, 0.000000001, sqrtx/SQRTPIO2, 0.0000001 );
    */

    return -2.0*SQRTPIO2*sqrtx*( ( 0.5 - c )*sin(x) - ( 0.5 - s )*cos(x) );
}



double fock_i( double x ) {

    double sqrtx, c, s;

    sqrtx = sqrt(x);

    fresnel( sqrtx/SQRTPIO2, &c, &s );

    /*
    s = quad_romberg( Sin, 0.000000001, sqrtx/SQRTPIO2, 0.0000001 );
    c = quad_romberg( Cos, 0.000000001, sqrtx/SQRTPIO2, 0.0000001 );
    */

    return 2.0*SQRTPIO2*sqrtx*( ( 0.5 - c )*cos(x) + ( 0.5 - s )*sin(x) );
}





dcomplex fock( double x ) {

    double sqrtx, c, s;
    dcomplex f;

    sqrtx = sqrt(x);

    fresnel( sqrt(x/PIO2), &c, &s );

    f.r = -2.0*SQRTPIO2*sqrtx*( ( 0.5 - c )*sin(x) - ( 0.5 - s )*cos(x) );
    f.i = 2.0*SQRTPIO2*sqrtx*( ( 0.5 - c )*cos(x) + ( 0.5 - s )*sin(x) );

    return f;
}



