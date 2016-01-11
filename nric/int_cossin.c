

#include <math.h>

#include "nric.h"


void cisi( double x, double *ci, double *si ) {

    const double EPS = 6.0e-10;
    const int MAXIT = 100;
    const int TRUE = 1;
    const double TMIN = 2.0;
    const double FPMIN = 1.0e-30;
    const dcomplex ONE = mk_dc( 1.0, 0.0);
    const dcomplex TWO = mk_dc( 2.0, 0.0);

    int i, k, odd;
    double a, err, fact, sign, sum, sumc, sums, t, term;
    dcomplex h, b, c, d, del;

    t = fabs(x);
    if ( t == 0.0 ) {
        *ci = -1.0/FPMIN;
        *si = 0.0;
        return;
    }
    if ( t > TMIN ) {
        /*
         *    Evaluate Ci and Si by Lentz's modified method of continued fracions.
         */
        b = mk_dc( 1.0, t );
        c = mk_dc( 1.0/FPMIN, 0.0 );
        d = h = div_dc( ONE, b );
        for ( i = 2; i <= MAXIT; ++i ) {
            a = -(i-1)*(i-1);
            b = add_dc( b, TWO );
            d = div_dc( ONE, add_dc( rmul_dc( a, d ), b ) );
            c = add_dc( b, div_dc( mk_dc( a, 0.0 ), c ) );
            del = mul_dc( c, d );
            h = mul_dc( h, del );
            if ( fabs( del.r - 1.0 ) + fabs( del.i ) < EPS ) break;
        }
        if ( i > MAXIT ) nrerror( "Continued fraction evaluation failed in cisi." );
        h = mul_dc( mk_dc( cos(t), -sin(t) ), h );
        *ci = -h.r;
        *si = PIO2 + h.i;
    } else {
        /*
         *    Evaluate Ci and Si by series simultaneously.
         */
        if ( t < sqrt(FPMIN) ) {
            /*
             *    Avoid underflow.
             */
            sumc = 0.0;
            sums = t;
        }  else {
            /*
             *    Evaluate S and C by series expansion.
             */
            sum = sums = sumc = 0.0;
            sign = fact = 1.0;
            odd = TRUE;
            for ( k = 1; k <= MAXIT; ++k ) {
                fact *= t/k;
                term = fact/k;
                sum += sign*term;
                err = term/fabs(sum);
                if ( odd ) {
                    sign = -sign;
                    sums = sum;
                    sum = sumc;
                } else {
                    sumc = sum;
                    sum = sums;
                }
                if ( err < EPS ) break;
                odd = !odd;
            }
            if ( k > MAXIT ) nrerror( "Series evaluation failed in cisi." );
        }
        *ci = sumc + log(t) + GAMMA;
        *si = sums;
    }
    if ( x < 0.0 ) *si = -*si;
}


