
#include <math.h>
#include <malloc.h>

#include "nric.h"


dcomplex *dcvector( int nl, int nh ) {

    dcomplex *v;

    v = ( dcomplex * ) malloc( ( unsigned int ) ( nh - nl + 1 )*sizeof( dcomplex ) );
    if ( !v ) nrerror( "Allocation failure in dcvector." );
    v -= nl;

    return  v;
}


void  free_dcvector( dcomplex *v, int nl, int nh ) {

    free( ( void * ) ( v + nl ) );
}


dcomplex **dcmatrix( int nrl, int nrh, int ncl, int nch ) {

    int i;
    dcomplex **m;

    /*
     *    Allocate pointers to rows.
     */
    m = ( dcomplex ** ) malloc( ( unsigned int ) ( nrh - nrl + 1 )*sizeof( dcomplex * ) );
    if ( !m ) nrerror( "Allocation failure 1 in dcmatrix." );
    m -= nrl;


    /*
     *    Allocate rows and set pointers to them.
     */
    for ( i = nrl; i <= nrh; ++i) {

        m[i] = ( dcomplex * ) malloc( ( unsigned int ) ( nch - ncl + 1 )*sizeof( dcomplex ) );
        if ( !m[i] ) nrerror( "Allocation failure 2 in dcmatrix." );
        m[i] -= ncl;
    }

    return m;
}


void  free_dcmatrix( dcomplex **m, int nrl, int nrh, int ncl, int nch ) {

    int i;

    for ( i = nrh; i >= nrl; --i ) free( ( void * ) ( m[i] + ncl ) );
    free( ( void * ) ( m + nrl ) );
}


dcomplex add_dc( dcomplex a, dcomplex b )  {
    dcomplex c;
    c.r = a.r + b.r;
    c.i = a.i + b.i;
    return c;
}



dcomplex sub_dc( dcomplex a, dcomplex b )  {
    dcomplex c;
    c.r = a.r - b.r;
    c.i = a.i - b.i;
    return c;
}



dcomplex mul_dc( dcomplex a, dcomplex b )  {
    dcomplex c;
    c.r = a.r*b.r - a.i*b.i;
    c.i = a.i*b.r + a.r*b.i;
    return c;
}



dcomplex div_dc( dcomplex a, dcomplex b )  {
    dcomplex c;
    double r, den;
    if ( fabs( b.r ) >= fabs( b.i ) ) {
        r = b.i/b.r;
        den = b.r + r*b.i;
        c.r = ( a.r + r*a.i )/den;
        c.i = ( a.i - r*a.r )/den;
    } else {
        r = b.r/b.i;
        den = b.i + r*b.r;
        c.r = ( r*a.r + a.i )/den;
        c.i = ( r*a.i - a.r )/den;
    }
    return c;
}



dcomplex mk_dc( double a, double b )  {
    dcomplex c;
    c.r = a;
    c.i = b;
    return c;
}



double abs_dc( dcomplex a )  {
    double x, y, r, temp;
    x = fabs( a.r );
    y = fabs( a.i );
    if ( x == 0.0 ) r = y;
    else if ( y == 0.0 ) r = x;
    else if ( x >= y ) {
        temp = y/x;
        r = x*sqrt( 1.0 + temp*temp );
    } else {
        temp = x/y;
        r = y*sqrt( 1.0 + temp*temp );
    }
    return r;
}


double phase_dc( dcomplex a )  {

    double p;
    p = atan2( a.i, a.r );
    return p;
}



dcomplex conj_dc( dcomplex a )  {
    dcomplex c;
    c.r = a.r;
    c.i = -a.i;
    return c;
}



dcomplex sqrt_dc( dcomplex a )  {
    dcomplex c;
    double x, y, w, r;
    if ( ( a.r == 0.0 ) && ( a.i == 0.0 ) )  {
        c.r = c.i = 0.0;
        return c;
    } else {
        x = fabs( a.r );
        y = fabs( a.i );
        if ( x >= y )  {
            r = y/x;
            w = sqrt(x)*sqrt( 0.5*( 1.0 + sqrt( 1.0 + r*r ) ) );
        } else {
            r = x/y;
            w = sqrt(y)*sqrt( 0.5*( r + sqrt( 1.0 + r*r ) ) );
        }
        if ( a.r >= 0.0 )  {
            c.r = w;
            c.i = a.i/(2.0*w);
        } else {
            c.i = (a.i >= 0.0) ? w : -w;
            c.r = a.i/(2.0*c.i);
        }
    }
    return c;
}



dcomplex iexp_dc( double a )  {
    dcomplex c;
    c.r = cos(a);
    c.i = sin(a);
    return c;
}



dcomplex exp_dc( dcomplex a )  {
    dcomplex c;
    c.r = exp( a.r )*cos( a.i );
    c.i = exp( a.r )*sin( a.i );
    return c;
}



dcomplex sin_dc( dcomplex a )  {
    dcomplex c;
    c.r = sin( a.r )*cosh( a.i );
    c.i = cos( a.i )*sinh( a.r );
    return c;
}


dcomplex cos_dc( dcomplex a )  {
    dcomplex c;
    c.r = cos( a.r )*cosh( a.i );
    c.i = -sin( a.i )*sinh( a.r );
    return c;
}



dcomplex tan_dc( dcomplex a )  {

    dcomplex s, c, t;
    double r, den;

    s.r = sin( a.r )*cosh( a.i );
    s.i = cos( a.i )*sinh( a.r );

    c.r = cos( a.r )*cosh( a.i );
    c.i = -sin( a.i )*sinh( a.r );

    if ( fabs( c.r ) >= fabs( c.i ) ) {
        r = c.i/c.r;
        den = c.r + r*c.i;
        t.r = ( s.r + r*s.i )/den;
        t.i = ( s.i - r*s.r )/den;
    } else {
        r = c.r/c.i;
        den = c.i + r*c.r;
        t.r = ( r*s.r + s.i )/den;
        t.i = ( r*s.i - s.r )/den;
    }
    return t;
}



dcomplex sinh_dc( dcomplex a )  {
    dcomplex c;
    c.r = sinh( a.r )*cos( a.i );
    c.i = cosh( a.r )*sin( a.i );
    return c;
}


dcomplex cosh_dc( dcomplex a )  {
    dcomplex c;
    c.r = cosh( a.r )*cos( a.i );
    c.i = sinh( a.r )*sin( a.i );
    return c;
}


dcomplex tanh_dc( dcomplex a )  {

    dcomplex s, c, t;
    double r, den;

    s.r = sinh( a.r )*cos( a.i );
    s.i = cosh( a.i )*sin( a.r );

    c.r = cosh( a.r )*cos( a.i );
    c.i = sinh( a.r )*sin( a.i );

    if ( fabs( c.r ) >= fabs( c.i ) ) {
        r = c.i/c.r;
        den = c.r + r*c.i;
        t.r = ( s.r + r*s.i )/den;
        t.i = ( s.i - r*s.r )/den;
    } else {
        r = c.r/c.i;
        den = c.i + r*c.r;
        t.r = ( r*s.r + s.i )/den;
        t.i = ( r*s.i - s.r )/den;
    }
    return t;
}


dcomplex log_dc( dcomplex a ) {
    double l = hypot( a.r, a.i );
    dcomplex c;
    c.r = log(l);
    c.i = atan2( a.i, a.r );
    return c;
}


dcomplex log10_dc( dcomplex a ) {
    double l = hypot( a.r, a.i );
    dcomplex c;
    c.r = log10(l);
    c.i = atan2( a.i, a.r )/log(10.);
    return c;
}



dcomplex rmul_dc( double a, dcomplex b )  {
    dcomplex c;
    c.r = a*b.r;
    c.i = a*b.i;
    return c;
}



dcomplex imul_dc( double a, dcomplex b )  {
    dcomplex c;
    c.r = a*b.i;
    c.i = a*b.r;
    return c;
}



