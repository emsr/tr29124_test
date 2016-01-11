
#include <math.h>


#include "nric.h"


double chebyshev_t0( double x ) {

    return 1.0;
}


double chebyshev_t1( double x ) {

    return x;
}


double chebyshev_t( int n, double x ) {

    double t0, t1, tm1, t, tp1;
    int j;

    t0 = 1.0;
    if ( n == 0 ) return t0;

    t1 = x;
    if ( n == 1 ) return t1;

    for ( tm1 = t0, t = t1, j = 1; j < n; ++j ) {
        tp1 = 2.0*x*t - tm1;
        tm1 = t;
        t = tp1;
    }
    return t;
}

