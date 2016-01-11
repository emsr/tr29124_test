

#include <math.h>

#include "nric.h"

/*
 *    This routine solves a quadratic equation:
 *    0 = a0 + a1.x + a2.x^2
 *    placing the (real) roots in r1 and r2
 *    or the real and imaginary parts of a complex root in r1 and r2 respectively.
 *    The number of unique real roots is returned.
 */
int quadratic( double a0, double a1, double a2, double *r1, double *r2 ) {

    double d, q, sign;

    d = a1*a1 - 4.0*a2*a0;

    if ( d < 0.0 ) {
        /*
         *    Set the first root to the real part of the true complex root
         *    and the second root to the imaginary part of the true complex root.
         */
        *r1 = -a1/(2.0*a2);
        *r2 = -dsign(sqrt(-d),a1)/(2.0*a2);
        return 0;
    } else if ( d == 0.0 ) {
        *r1 = *r2 = -a1/(2.0*a2);
        return 1;
    } else {
        q = -0.5*(a1 + dsign(sqrt(d),a1));
        *r1 = q/a2;
        *r2 = a0/q;
        return 2;
    }
}

/*
 *    This routine solves a quadratic equation:
 *    0 = a0 + a1.x + a2.x^2
 *    placing the (complex) roots in r1 and r2.
 */
void quad_complex( dcomplex a0, dcomplex a1, dcomplex a2, dcomplex *r1, dcomplex *r2 ) {

    double sign;
    dcomplex d, q;

    d = sub_dc( mul_dc( a1, a1 ), rmul_dc( 4.0, mul_dc( a2, a0) ) );

    sign = mul_dc( conj_dc(a1), d ).r > 0.0 ? +1.0 : -1.0;
    q = rmul_dc( -0.5, add_dc( a1, rmul_dc( sign, sqrt_dc( d ) ) ) );
    *r1 = div_dc( q, a2);
    *r2 = div_dc( a0, q);

    return;
}

