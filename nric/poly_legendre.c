

#include <math.h>


#include "nric.h"



double legendre_p0( double x ) {

    if ( ( x < -1.000001 ) || ( x > +1.000001 ) ) nrerror( "Bad argument x in legendre_p0." );

    return 1.0;
}



double legendre_p1( double x ) {

    if ( ( x < -1.00000001 ) || ( x > +1.00000001 ) ) nrerror( "Bad argument x in legendre_p1." );

    return x;
}



double legendre_p( int l, double x ) {

    int j;
    double pm2, pm1, p;

    if ( l < 0 ) nrerror( "Bad order in legendre_p." );
    if ( ( x < -1.00000001 ) || ( x > +1.00000001 ) ) nrerror( "Bad argument x in legendre_p." );

    pm2 = 1.0;
    if ( l == 0 ) return pm2;
    pm1 = x;
    if ( l == 1 ) return pm1;
    for ( j = 2; j <= l; ++j ) {
        p = ( ( 2.0*l - 1.0 )*x*pm1 - ( 1.0*l - 1.0 )*pm2 )/l;
        pm2 = pm1;
        pm1 = p;
    }
    return p;
}


double legendre_q0( double x ) {

    if ( ( x < -0.999999999 ) || ( x > +0.999999999 ) ) nrerror( "Bad argument x in legendre_p0." );

    return 0.5*log((1.0+x)/(1.0-x));
}


double legendre_q1( double x ) {

    if ( ( x < -0.999999999 ) || ( x > +0.999999999 ) ) nrerror( "Bad argument x in legendre_p0." );

    return 0.5*x*log((1.0+x)/(1.0-x)) - 1.0;
}



double legendre_q( int l, double x ) {

    int j;
    double qm2, qm1, q;

    if ( l < 0 ) nrerror( "Bad order in legendre_q." );
    if ( ( x < -0.999999999 ) || ( x > +0.999999999 ) ) nrerror( "Bad argument x in legendre_q." );

    qm2 = 0.5*log((1.0+x)/(1.0-x));
    if ( l == 0 ) return qm2;
    qm1 = 0.5*x*log((1.0+x)/(1.0-x)) - 1.0;
    if ( l == 1 ) return qm1;
    for ( j = 2; j <= l; ++j ) {
        q = ( ( 2.0*l - 1.0 )*x*qm1 - ( 1.0*l - 1.0 )*qm2 )/l;
        qm2 = qm1;
        qm1 = q;
    }
    return q;
}


