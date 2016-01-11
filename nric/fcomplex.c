
#include <malloc.h>
#include <math.h>

#include "nric.h"


fcomplex *fcvector( int nl, int nh ) {

    fcomplex *v;

    v = ( fcomplex * ) malloc( ( unsigned int ) ( nh - nl + 1 )*sizeof( fcomplex ) );
    if ( !v ) nrerror( "Allocation failure in fcvector." );
    v -= nl;

    return  v;
}


void  free_fcvector( fcomplex *v, int nl, int nh ) {

    free( ( void * ) ( v + nl ) );
}


fcomplex **fcmatrix( int nrl, int nrh, int ncl, int nch ) {

    int i;
    fcomplex **m;

    /*
     *    Allocate pointers to rows.
     */
    m = ( fcomplex ** ) malloc( ( unsigned int ) ( nrh - nrl + 1 )*sizeof( fcomplex * ) );
    if ( !m ) nrerror( "Allocation failure 1 in fcmatrix." );
    m -= nrl;


    /*
     *    Allocate rows and set pointers to them.
     */
    for ( i = nrl; i <= nrh; ++i) {

        m[i] = ( fcomplex * ) malloc( ( unsigned int ) ( nch - ncl + 1 )*sizeof( fcomplex ) );
        if ( !m[i] ) nrerror( "Allocation failure 2 in fcmatrix." );
        m[i] -= ncl;
    }

    return m;
}


void  free_fcmatrix( fcomplex **m, int nrl, int nrh, int ncl, int nch ) {

    int i;

    for ( i = nrh; i >= nrl; --i ) free( ( void * ) ( m[i] + ncl ) );
    free( ( void * ) ( m + nrl ) );
}


fcomplex add_fc( fcomplex a, fcomplex b )  {
    fcomplex c;
    c.r = a.r + b.r;
    c.i = a.i + b.i;
    return c;
}


fcomplex sub_fc( fcomplex a, fcomplex b )  {
    fcomplex c;
    c.r = a.r - b.r;
    c.i = a.i - b.i;
    return c;
}



fcomplex mul_fc( fcomplex a, fcomplex b )  {
    fcomplex c;
    c.r = a.r*b.r - a.i*b.i;
    c.i = a.i*b.r + a.r*b.i;
    return c;
}



fcomplex div_fc( fcomplex a, fcomplex b )  {
    fcomplex c;
    float r, den;
    if ( fabsf( b.r ) >= fabsf( b.i ) ) {
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



fcomplex mk_fc( float a, float b )  {
    fcomplex c;
    c.r = a;
    c.i = b;
    return c;
}



float abs_fc( fcomplex a )  {
    float x, y, r, temp;
    x = fabsf( a.r );
    y = fabsf( a.i );
    if ( x == 0.0 ) r = y;
    else if ( y == 0.0 ) r = x;
    else if ( x >= y ) {
        temp = y/x;
        r = x*sqrtf( 1.0 + temp*temp );
    } else {
        temp = x/y;
        r = y*sqrtf( 1.0 + temp*temp );
    }
    return r;
}


float phase_fc( fcomplex a )  {

    float p;
    p = atan2f( a.i, a.r );
    return p;
}



fcomplex conj_fc( fcomplex a )  {
    fcomplex c;
    c.r = a.r;
    c.i = -a.i;
    return c;
}



fcomplex sqrt_fc( fcomplex a )  {
    fcomplex c;
    float x, y, w, r;
    if ( ( a.r == 0.0 ) && ( a.i == 0.0 ) )  {
        c.r = c.i = 0.0;
        return c;
    } else {
        x = fabsf( a.r );
        y = fabsf( a.i );
        if ( x >= y )  {
            r = y/x;
            w = sqrtf(x)*sqrtf( 0.5*( 1.0 + sqrtf( 1.0 + r*r ) ) );
        } else {
            r = x/y;
            w = sqrtf(y)*sqrtf( 0.5*( r + sqrtf( 1.0 + r*r ) ) );
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



fcomplex iexp_fc( float a )  {
    fcomplex c;
    c.r = cosf(a);
    c.i = sinf(a);
    return c;
}



fcomplex exp_fc( fcomplex a )  {
    fcomplex c;
    c.r = expf( a.r )*cosf( a.i );
    c.i = expf( a.r )*sinf( a.i );
    return c;
}



fcomplex sin_fc( fcomplex a )  {
    fcomplex c;
    c.r = sinf( a.r )*coshf( a.i );
    c.i = cosf( a.i )*sinhf( a.r );
    return c;
}


fcomplex cos_fc( fcomplex a )  {
    fcomplex c;
    c.r = cosf( a.r )*coshf( a.i );
    c.i = -sinf( a.i )*sinhf( a.r );
    return c;
}



fcomplex tan_fc( fcomplex a )  {

    fcomplex s, c, t;
    float r, den;

    s.r = sinf( a.r )*coshf( a.i );
    s.i = cosf( a.i )*sinhf( a.r );

    c.r = cosf( a.r )*coshf( a.i );
    c.i = -sinf( a.i )*sinhf( a.r );

    if ( fabsf( c.r ) >= fabsf( c.i ) ) {
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



fcomplex sinh_fc( fcomplex a )  {
    fcomplex c;
    c.r = sinhf( a.r )*cosf( a.i );
    c.i = coshf( a.r )*sinf( a.i );
    return c;
}


fcomplex cosh_fc( fcomplex a )  {
    fcomplex c;
    c.r = coshf( a.r )*cosf( a.i );
    c.i = sinhf( a.r )*sinf( a.i );
    return c;
}


fcomplex tanh_fc( fcomplex a )  {

    fcomplex s, c, t;
    float r, den;

    s.r = sinhf( a.r )*cosf( a.i );
    s.i = coshf( a.i )*sinf( a.r );

    c.r = coshf( a.r )*cosf( a.i );
    c.i = sinhf( a.r )*sinf( a.i );

    if ( fabsf( c.r ) >= fabsf( c.i ) ) {
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


fcomplex log_fc( fcomplex a ) {
    float l = hypot( a.r, a.i );
    fcomplex c;
    c.r = log(l);
    c.i = atan2( a.i, a.r );
    return c;
}


fcomplex log10_fc( fcomplex a ) {
    float l = hypot( a.r, a.i );
    fcomplex c;
    c.r = log10(l);
    c.i = atan2( a.i, a.r )/log(10.);
    return c;
}



fcomplex rmul_fc( float a, fcomplex b )  {
    fcomplex c;
    c.r = a*b.r;
    c.i = a*b.i;
    return c;
}



fcomplex imul_fc( float a, fcomplex b )  {
    fcomplex c;
    c.r = a*b.i;
    c.i = a*b.r;
    return c;
}



