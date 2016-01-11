
#include <math.h>

#include "nric.h"

#define NMAX 10
#define XMIN 0.4
#define H 0.2

double dawson( double x ) {

    int i, n0;
    double d1, d2, e1, e2, sum, x2, xp, xx, ans;
    static double c[NMAX+1];
    static int init = 0;

    if ( init == 0 ) {
        init = 1;
        for ( i = 1; i <= NMAX; ++i ) {c[i] = exp( -dsqr( ( 2.0*i - 1.0 )*H ) );}
    }
    if ( fabs(x) < XMIN ) {
        /*
         *    Use series expansion.
         */
        x2 = x*x;
        ans = x*( 1.0 - (2.0/3.0)*x2*( 1.0 - (2.0/5.0)*x2*( 1.0 - (2.0/7.0)*x2 ) ) );
    } else {
        /*
         *    Use sampling theorem representation.
         */
        xx = fabs(x);
        n0 = 2*( (int)(0.5 + 0.5*xx/H) );
        xp = xx - n0*H;
        e1 = exp(2.0*xp*H);
        e2 = e1*e1;
        d1 = (double) n0 + 1.0;
        d2 = d1 - 2.0;
        sum = 0.0;
        for ( i = 1; i <= NMAX; ++i ) {
            sum += c[i]*( e1/d1 + 1.0/(d2*e1) );
            d1 += 2.0;
            d2 -= 2.0;
            e1 *= e2;
        }
        ans = dsign( exp( -xp*xp ), x )*sum/SQRTPI;
    }
    return ans;
}



