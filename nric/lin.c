
#include <math.h>
#include <stdio.h>

#include "nric.h"

int main() {

    double *v1, *v2, *v;
    double **m1, **m2, **m, **ms, **ma, **mt, **mo;
    double **mtemp;
    double *vtemp;
    double xtemp;
    double trace;

    v1 = dvector( 1, 2 );
    v2 = dvector( 1, 2 );
    v = dvector( 1, 2 );
    vtemp = dvector( 1, 2 );

    m1 = dmatrix( 1, 2, 1, 2);
    m2 = dmatrix( 1, 2, 1, 2);
    m = dmatrix( 1, 2, 1, 2);
    ms = dmatrix( 1, 2, 1, 2);
    ma = dmatrix( 1, 2, 1, 2);
    mt = dmatrix( 1, 2, 1, 2);
    mo = dmatrix( 1, 2, 1, 2);
    mtemp = dmatrix( 1, 2, 1, 2);

    m1[1][1] = 1.0;
    m1[1][2] = 2.0;
    m1[2][1] = 3.0;
    m1[2][2] = 4.0;

    v1[1] = 5.0;
    v1[2] = 6.0;

    printf( "Hello 1\n" );

    dcopyvec( v1, 1, 2, v2 );
    dcopymat( m1, 1, 2, 1, 2, m2 );

    printf( "Hello 2\n" );

    dvecvec( v1, v2, 1, 2, &xtemp );
    printf( " v1.v2 = x = %4.1f\n", xtemp );

    dmatvec( m1, v2, 1, 2, 1, 2, vtemp );
    printf( " m1.v2 = v = %4.1f, %4.1f\n", vtemp[1], vtemp[2] );

    dvecmat( v1, m2, 1, 2, 1, 2, vtemp );
    printf( " v1.m2 = v = %4.1f, %4.1f\n", vtemp[1], vtemp[2] );

    dmatmat( m1, m2, 1, 2, 1, 2, 1, 2, mtemp );
    printf( " m1.m2 = x = %4.1f, %4.1f, %4.1f, %4.1f\n", mtemp[1][1], mtemp[1][2], mtemp[2][1], mtemp[2][2] );

    daddvec( v1, v2, 1, 2, v );
    printf( " v1 + v2 = v = %4.1f, %4.1f\n", v[1], v[2] );

    daddmat( m1, m2, 1, 2, 1, 2, m );
    printf( " m1 + m2 = m = %4.1f, %4.1f, %4.1f, %4.1f\n", m[1][1], m[1][2], m[2][1], m[2][2] );

    dsubvec( v1, v2, 1, 2, v );
    printf( " v1 - v2 = v = %4.1f, %4.1f\n", v[1], v[2] );

    dsubmat( m1, m2, 1, 2, 1, 2, m );
    printf( " m1 - m2 = m = %4.1f, %4.1f, %4.1f, %4.1f\n", m[1][1], m[1][2], m[2][1], m[2][2] );

    dmulvec( 2.0, v2, 1, 2, v );
    printf( " 2*v2 = v = %4.1f, %4.1f\n", v[1], v[2] );

    dmulmat( 2.0, m2, 1, 2, 1, 2, m );
    printf( " 2*m2 = m = %4.1f, %4.1f, %4.1f, %4.1f\n", m[1][1], m[1][2], m[2][1], m[2][2] );

    dtracemat( m1, 1, 2, &trace );
    printf( " trace = %5.1f\n", trace );

    dtransmat( m1, 1, 2, 1, 2, mt );
    printf( " mt = %4.1f, %4.1f, %4.1f, %4.1f\n", mt[1][1], mt[1][2], mt[2][1], mt[2][2] );

    dsymmat( m1, 1, 2, ms, ma );
    printf( " ms = %4.1f, %4.1f, %4.1f, %4.1f\n", ms[1][1], ms[1][2], ms[2][1], ms[2][2] );
    printf( " ma = %4.1f, %4.1f, %4.1f, %4.1f\n", ma[1][1], ma[1][2], ma[2][1], ma[2][2] );

    dgramm_schmidt( m, 1, 2, mo );
    printf( " m  = %4.1f, %4.1f, %4.1f, %4.1f\n", m[1][1], m[1][2], m[2][1], m[2][2] );
    printf( " mo = %4.1f, %4.1f, %4.1f, %4.1f\n", mo[1][1], mo[1][2], mo[2][1], mo[2][2] );
}


