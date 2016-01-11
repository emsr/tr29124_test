
#include <math.h>
#include <malloc.h>

#include "nric.h"


/************************************************************************************
 ************************************************************************************
 **
 **    MIXED routines.
 **
 ************************************************************************************
 ************************************************************************************/

fmatrix32 cart_vec3vec2_fv( fvector3 U, fvector2 V ) {
    fmatrix32 m;

    m.x.x = U.x*V.x;
    m.x.y = U.x*V.y;
    m.y.x = U.y*V.x;
    m.y.y = U.y*V.y;
    m.z.x = U.z*V.x;
    m.z.y = U.z*V.y;

    return m;
}


fmatrix23 cart_vec2vec3_fv( fvector2 U, fvector3 V ) {
    fmatrix23 m;

    m.x.x = U.x*V.x;
    m.x.y = U.x*V.y;
    m.x.z = U.x*V.z;
    m.y.x = U.y*V.x;
    m.y.y = U.y*V.y;
    m.y.z = U.y*V.z;

    return m;
}

fvector3 mat32vec2_fv( fmatrix32 A, fvector2 V ) {
    fvector3 v;

    v.x = A.x.x*V.x + A.x.y*V.y;
    v.y = A.y.x*V.x + A.y.y*V.y;
    v.z = A.z.x*V.x + A.z.y*V.y;

    return v;
}

fvector2 mat23vec3_fv( fmatrix23 A, fvector3 V ) {
    fvector2 v;

    v.x = A.x.x*V.x + A.x.y*V.y + A.x.z*V.z;
    v.y = A.y.x*V.x + A.y.y*V.y + A.y.z*V.z;

    return v;
}

fmatrix32 trans_mat23_fv( fmatrix23 A ) {
    fmatrix32 C;

    C.x.x = A.x.x;
    C.x.y = A.y.x;
    C.y.x = A.x.y;
    C.y.y = A.y.y;
    C.z.x = A.x.z;
    C.z.y = A.y.z;

    return C;
}

fmatrix23 trans_mat32_fv( fmatrix32 A ) {
    fmatrix23 C;

    C.x.x = A.x.x;
    C.x.y = A.y.x;
    C.x.z = A.z.x;
    C.y.x = A.x.y;
    C.y.y = A.y.y;
    C.y.z = A.z.y;

    return C;
}


fmatrix3 mat32mat23_fv( fmatrix32 A, fmatrix23 B ) {
    fmatrix3 C;

    C.x.x = A.x.x*B.x.x + A.x.y*B.y.x;
    C.x.y = A.x.x*B.x.y + A.x.y*B.y.y;
    C.x.y = A.x.x*B.x.z + A.x.y*B.y.z;
    C.y.x = A.y.x*B.x.x + A.y.y*B.y.x;
    C.y.y = A.y.x*B.x.y + A.y.y*B.y.y;
    C.y.z = A.y.x*B.x.z + A.y.y*B.y.z;
    C.z.x = A.z.x*B.x.x + A.z.y*B.y.x;
    C.z.y = A.z.x*B.x.y + A.z.y*B.y.y;
    C.z.z = A.z.x*B.x.z + A.z.y*B.y.z;

    return C;
}


fmatrix32 mat33mat32_fv( fmatrix3 A, fmatrix32 B ) {
    fmatrix32 C;

    C.x.x = A.x.x*B.x.x + A.x.y*B.y.x + A.x.z*B.z.x;
    C.x.y = A.x.x*B.x.y + A.x.y*B.y.y + A.x.z*B.z.y;
    C.y.x = A.y.x*B.x.x + A.y.y*B.y.x + A.y.z*B.z.x;
    C.y.y = A.y.x*B.x.y + A.y.y*B.y.y + A.y.z*B.z.y;
    C.z.x = A.z.x*B.x.x + A.z.y*B.y.x + A.z.z*B.z.x;
    C.z.y = A.z.x*B.x.y + A.z.y*B.y.y + A.z.z*B.z.y;

    return C;
}


fmatrix32 mat32mat22_fv( fmatrix32 A, fmatrix2 B ) {
    fmatrix32 C;

    C.x.x = A.x.x*B.x.x + A.x.y*B.y.x;
    C.x.y = A.x.x*B.x.y + A.x.y*B.y.y;
    C.y.x = A.y.x*B.x.x + A.y.y*B.y.x;
    C.y.y = A.y.x*B.x.y + A.y.y*B.y.y;
    C.z.x = A.z.x*B.x.x + A.z.y*B.y.x;
    C.z.y = A.z.x*B.x.y + A.z.y*B.y.y;

    return C;
}

fmatrix2 mat23mat32_fv( fmatrix23 A, fmatrix32 B ) {
    fmatrix2 C;

    C.x.x = A.x.x*B.x.x + A.x.y*B.y.x + A.x.z*B.z.x;
    C.x.y = A.x.x*B.x.y + A.x.y*B.y.y + A.x.z*B.z.y;
    C.y.x = A.y.x*B.x.x + A.y.y*B.y.x + A.y.z*B.z.x;
    C.y.y = A.y.x*B.x.y + A.y.y*B.y.y + A.y.z*B.z.y;

    return C;
}


fmatrix23 mat22mat23_fv( fmatrix2 A, fmatrix23 B ) {
    fmatrix23 C;

    C.x.x = A.x.x*B.x.x + A.x.y*B.y.x;
    C.x.y = A.x.x*B.x.y + A.x.y*B.y.y;
    C.x.z = A.x.x*B.x.z + A.x.y*B.y.z;
    C.y.x = A.y.x*B.x.x + A.y.y*B.y.x;
    C.y.y = A.y.x*B.x.y + A.y.y*B.y.y;
    C.y.z = A.y.x*B.x.z + A.y.y*B.y.z;

    return C;
}


fmatrix23 mat23mat33_fv( fmatrix23 A, fmatrix3 B ) {
    fmatrix23 C;

    C.x.x = A.x.x*B.x.x + A.x.y*B.y.x + A.x.z*B.z.x;
    C.x.y = A.x.x*B.x.y + A.x.y*B.y.y + A.x.z*B.z.y;
    C.x.z = A.x.x*B.x.z + A.x.y*B.y.z + A.x.z*B.z.z;
    C.y.x = A.y.x*B.x.x + A.y.y*B.y.x + A.y.z*B.z.x;
    C.y.y = A.y.x*B.x.y + A.y.y*B.y.y + A.y.z*B.z.y;
    C.y.z = A.y.x*B.x.z + A.y.y*B.y.z + A.y.z*B.z.z;

    return C;
}




/************************************************************************************
 ************************************************************************************
 **
 **    VECTOR3 routines.
 **
 ************************************************************************************
 ************************************************************************************/


fvector3 *fv3vector( int nl, int nh ) {

    fvector3 *v;

    v = ( fvector3 * ) malloc( ( unsigned int ) ( nh - nl + 1 )*sizeof( fvector3 ) );
    if ( !v ) nrerror( "Allocation failure in fv3vector." );
    v -= nl;

    return  v;
}

void  free_fv3vector( fvector3 *v, int nl, int nh ) {
    free( ( void * ) ( v + nl ) );
}


fvector2 *fv2vector( int nl, int nh ) {

    fvector2 *v;

    v = ( fvector2 * ) malloc( ( unsigned int ) ( nh - nl + 1 )*sizeof( fvector2 ) );
    if ( !v ) nrerror( "Allocation failure in fv2vector." );
    v -= nl;

    return  v;
}

void  free_fv2vector( fvector2 *v, int nl, int nh ) {
    free( ( void * ) ( v + nl ) );
}



float *convert_fvector_fv3( fvector3 a ) {

    float *v;
    v = fvector( 1, 3 );
    v[1] = a.x;
    v[2] = a.y;
    v[3] = a.z;
    return v;
}

void free_convert_fvector_fv3( float *v ) {

    free_fvector( v, 1, 3 );
}



float **convert_fmatrix_fv3( fmatrix3 a ) {

    float **m;
    m = fmatrix( 1, 3, 1, 3 );
    m[1][1] = a.x.x;
    m[1][2] = a.x.y;
    m[1][3] = a.x.z;
    m[2][1] = a.y.x;
    m[2][2] = a.y.y;
    m[2][3] = a.y.z;
    m[3][1] = a.z.x;
    m[3][2] = a.z.y;
    m[3][3] = a.z.z;
    return m;
}

void free_convert_fmatrix_fv3( float **m ) {

    free_fmatrix( m, 1, 3, 1, 3 );
}


fvector3 mk_fvector_fv3( float *v ) {

    fvector3 a;
    a.x = v[1];
    a.y = v[2];
    a.z = v[3];
    return a;
}


fmatrix3 mk_fmatrix_fv3( float **m ) {

    fmatrix3 a;
    a.x.x = m[1][1];
    a.x.y = m[1][2];
    a.x.z = m[1][3];
    a.y.x = m[2][1];
    a.y.y = m[2][2];
    a.y.z = m[2][3];
    a.z.x = m[3][1];
    a.z.y = m[3][2];
    a.z.z = m[3][3];
    return a;
}


fvector3 mk_fv3( float x, float y, float z ) {

    fvector3 c;
    c.x = x;
    c.y = y;
    c.z = z;
    return c;
}



fvector3 add_fv3( fvector3 a, fvector3 b ) {

    fvector3 c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    return c;
}




fvector3 sub_fv3( fvector3 a, fvector3 b ) {

    fvector3 c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    return c;
}




float mod_fv3( fvector3 a ) {

    float max, myx, rx, ry, rz;

    max = ((rz = fabsf(a.z)) >= (myx = ((ry = fabsf(a.y)) >= (rx = fabsf(a.x)) ? ry : rx)) ? rz : myx );

    if ( max == 0.0 ) return 0.0;
    else {
        rx /= max;
        ry /= max;
        rz /= max;
        return max*sqrtf( rx*rx + ry*ry + rz*rz );
    }
}





fvector3 norm_fv3( fvector3 a ) {

    float max, myx, rx, ry, rz;

    max = ((rz = fabsf(a.z)) >= (myx = ((ry = fabsf(a.y)) >= (rx = fabsf(a.x)) ? ry : rx)) ? rz : myx );

    if ( max == 0.0 ) {
        a.x = a.y = a.z = 0.0;
    } else {
        rx /= max;
        ry /= max;
        rz /= max;
        max *= sqrtf( rx*rx + ry*ry + rz*rz );

        a.x /= max;
        a.y /= max;
        a.z /= max;
    }
    return a;
}




float dot_fv3( fvector3 a, fvector3 b ) {

    float c;
    c = a.x*b.x + a.y*b.y + a.z*b.z;
    return c;
}




float az_fv3( fvector3 a ) {

    float c;
    c = atan2f( a.y, a.x );
    return c;
}




float el_fv3( fvector3 a ) {

    float c;
    c = atan2f( a.z, sqrtf( a.x*a.x + a.y*a.y ) );
    return c;
}




fvector3 cross_fv3( fvector3 a, fvector3 b ) {

    fvector3 c;
    c.x = a.y*b.z - a.z*b.y;
    c.y = a.z*b.x - a.x*b.z;
    c.z = a.x*b.y - a.y*b.x;
    return c;
}




fvector3 rmul_fv3( float a, fvector3 b ) {

    fvector3 c;
    c.x = a*b.x;
    c.y = a*b.y;
    c.z = a*b.z;
    return c;
}




fmatrix3 cart_fv3( fvector3 a, fvector3 b ) {

    fmatrix3 c;
    c.x.x = a.x*b.x;
    c.x.y = a.x*b.y;
    c.x.z = a.x*b.z;
    c.y.x = a.y*b.x;
    c.y.y = a.y*b.y;
    c.y.z = a.y*b.z;
    c.z.x = a.z*b.x;
    c.z.y = a.z*b.y;
    c.z.z = a.z*b.z;
    return c;
}




fmatrix3 trans_fv3( fmatrix3 A ) {

    fmatrix3 a;

    a.x.x = A.x.x;
    a.x.y = A.y.x;
    a.x.z = A.z.x;
    a.y.x = A.x.y;
    a.y.y = A.y.y;
    a.y.z = A.z.y;
    a.z.x = A.x.z;
    a.z.y = A.y.z;
    a.z.z = A.z.z;

    return a;
}



/******************************************************************
 *
 *  VECTOR2 routines.
 *
 ******************************************************************/



float *convert_fvector_fv2( fvector2 a ) {

    float *v;
    v = fvector( 1, 2 );
    v[1] = a.x;
    v[2] = a.y;
    return v;
}

void free_convert_fvector_fv2( float *v ) {

    free_fvector( v, 1, 2 );
}


float **convert_fmatrix_fv2( fmatrix2 a ) {

    float **m;
    m = fmatrix( 1, 2, 1, 2 );
    m[1][1] = a.x.x;
    m[1][2] = a.x.y;
    m[2][1] = a.y.x;
    m[2][2] = a.y.y;
    return m;
}

void free_convert_fmatrix_fv2( float **m ) {

    free_fmatrix( m, 1, 2, 1, 2 );
}


fvector2 mk_fvector_fv2( float *v ) {

    fvector2 a;
    a.x = v[1];
    a.y = v[2];
    return a;
}


fmatrix2 mk_fmatrix_fv2( float **m ) {

    fmatrix2 a;
    a.x.x = m[1][1];
    a.x.y = m[1][2];
    a.y.x = m[2][1];
    a.y.y = m[2][2];
    return a;
}


fvector2 mk_fv2( float x, float y ) {

    fvector2 c;
    c.x = x;
    c.y = y;
    return c;
}



fvector2 add_fv2( fvector2 a, fvector2 b ) {

    fvector2 c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    return c;
}




fvector2 sub_fv2( fvector2 a, fvector2 b ) {

    fvector2 c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    return c;
}




float mod_fv2( fvector2 a ) {

    float max, rx, ry;

    max = ((ry = fabsf(a.y)) >= (rx = fabsf(a.x)) ? ry : rx);

    if ( max == 0.0 ) return 0.0;
    else {
        rx /= max;
        ry /= max;
        return max*sqrtf( rx*rx + ry*ry );
    }
}





fvector2 norm_fv2( fvector2 a ) {

    float max, rx, ry;

    max = ((ry = fabsf(a.y)) >= (rx = fabsf(a.x)) ? ry : rx);

    if ( max == 0.0 ) {
        a.x = a.y = 0.0;
    } else {
        rx /= max;
        ry /= max;
        max *= sqrtf( rx*rx + ry*ry );

        a.x /= max;
        a.y /= max;
    }
    return a;
}




float dot_fv2( fvector2 a, fvector2 b ) {

    float c;
    c = a.x*b.x + a.y*b.y ;
    return c;
}




float az_fv2( fvector2 a ) {

    float c;
    c = atan2f( a.y, a.x );
    return c;
}



fvector2 rmul_fv2( float a, fvector2 b ) {

    fvector2 c;
    c.x = a*b.x;
    c.y = a*b.y;
    return c;
}




fmatrix2 cart_fv2( fvector2 a, fvector2 b ) {

    fmatrix2 c;
    c.x.x = a.x*b.x;
    c.x.y = a.x*b.y;
    c.y.x = a.y*b.x;
    c.y.y = a.y*b.y;
    return c;
}




fmatrix2 trans_fv2( fmatrix2 A ) {

    fmatrix2 a;

    a.x.x = A.x.x;
    a.x.y = A.y.x;
    a.y.x = A.x.y;
    a.y.y = A.y.y;

    return a;
}



float trace_fv3( fmatrix3 a ) { return a.x.x + a.y.y + a.z.z; }

float det_fv3( fmatrix3 A ) {

     return A.x.x*(A.y.y*A.z.z - A.y.z*A.z.y)
          + A.x.y*(A.y.z*A.z.x - A.y.x*A.z.z)
          + A.x.z*(A.y.x*A.z.y - A.y.y*A.z.x);
}

fmatrix3 inv_fv3( fmatrix3 A ) {

    float det;
    fmatrix3 a;

    det = A.x.x*(A.y.y*A.z.z - A.y.z*A.z.y)
        + A.x.y*(A.y.z*A.z.x - A.y.x*A.z.z)
        + A.x.z*(A.y.x*A.z.y - A.y.y*A.z.x);

    a.x.x = (A.y.y*A.z.z - A.z.y*A.y.z)/det;
    a.y.x = (A.y.z*A.z.x - A.z.z*A.y.x)/det;
    a.z.x = (A.y.x*A.z.y - A.z.x*A.y.y)/det;
    a.x.y = (A.z.y*A.x.z - A.x.y*A.z.z)/det;
    a.y.y = (A.z.z*A.x.x - A.x.z*A.z.x)/det;
    a.z.y = (A.z.x*A.x.y - A.x.x*A.z.y)/det;
    a.x.z = (A.x.y*A.y.z - A.y.y*A.x.z)/det;
    a.y.z = (A.x.z*A.y.x - A.y.z*A.x.x)/det;
    a.z.z = (A.x.x*A.y.y - A.y.x*A.x.y)/det;

    return a;
}



float trace_fv2( fmatrix2 A ) { return A.x.x + A.y.y; }


float det_fv2( fmatrix2 A ) { return A.x.x*A.y.y - A.x.y*A.y.x; }


fmatrix2 inv_fv2( fmatrix2 A ) {

    float det;
    fmatrix2 a;

    det = A.x.x*A.y.y - A.x.y*A.y.x;

    a.x.x = A.y.y/det;
    a.y.x = -A.y.x/det;
    a.x.y = A.x.y/det;
    a.y.y = -A.x.x/det;

    return a;
}



fvector3 matvec_fv3( fmatrix3 A, fvector3 X ) {

    fvector3 V;

    V.x = A.x.x*X.x + A.x.y*X.y + A.x.z*X.z;
    V.y = A.y.x*X.x + A.y.y*X.y + A.y.z*X.z;
    V.z = A.z.x*X.x + A.z.y*X.y + A.z.z*X.z;

    return V;
}



fvector3 vecmat_fv3( fvector3 Y, fmatrix3 A ) {

    fvector3 V;

    V.x = Y.x*A.x.x + Y.y*A.y.x + Y.z*A.z.x;
    V.y = Y.x*A.x.y + Y.y*A.y.y + Y.z*A.z.y;
    V.z = Y.x*A.x.z + Y.y*A.y.z + Y.z*A.z.z;

    return V;
}



float vecmatvec_fv3( fvector3 Y, fmatrix3 A, fvector3 X ) {

    return (Y.x*A.x.x + Y.y*A.y.x + Y.z*A.z.x)*X.x
         + (Y.x*A.x.y + Y.y*A.y.y + Y.z*A.z.y)*X.y
         + (Y.x*A.x.z + Y.y*A.y.z + Y.z*A.z.z)*X.z;
}








fvector2 matvec_fv2( fmatrix2 A, fvector2 X ) {

    fvector2 V;

    V.x = A.x.x*X.x + A.x.y*X.y;
    V.y = A.y.x*X.x + A.y.y*X.y;

    return V;
}



fvector2 vecmat_fv2( fvector2 Y, fmatrix2 A ) {

    fvector2 V;

    V.x = Y.x*A.x.x + Y.y*A.y.x;
    V.y = Y.x*A.x.y + Y.y*A.y.y;

    return V;
}



float vecmatvec_fv2( fvector2 Y, fmatrix2 A, fvector2 X ) {

    return (Y.x*A.x.x + Y.y*A.y.x)*X.x
         + (Y.x*A.x.y + Y.y*A.y.y)*X.y;
}









fmatrix3 matmat_fv3( fmatrix3 A, fmatrix3 B ) {

    fmatrix3 C;

    C.x.x = A.x.x*B.x.x + A.x.y*B.y.x + A.x.z*B.z.x;
    C.x.y = A.x.x*B.x.y + A.x.y*B.y.y + A.x.z*B.z.y;
    C.x.z = A.x.x*B.x.z + A.x.y*B.y.z + A.x.z*B.z.z;
    C.y.x = A.y.x*B.x.x + A.y.y*B.y.x + A.y.z*B.z.x;
    C.y.y = A.y.x*B.x.y + A.y.y*B.y.y + A.y.z*B.z.y;
    C.y.z = A.y.x*B.x.z + A.y.y*B.y.z + A.y.z*B.z.z;
    C.z.x = A.z.x*B.x.x + A.z.y*B.y.x + A.z.z*B.z.x;
    C.z.y = A.z.x*B.x.y + A.z.y*B.y.y + A.z.z*B.z.y;
    C.z.z = A.z.x*B.x.z + A.z.y*B.y.z + A.z.z*B.z.z;

    return C;
}



fmatrix2 matmat_fv2( fmatrix2 A, fmatrix2 B ) {

    fmatrix2 C;

    C.x.x = A.x.x*B.x.x + A.x.y*B.y.x;
    C.x.y = A.x.x*B.x.y + A.x.y*B.y.y;
    C.y.x = A.y.x*B.x.x + A.y.y*B.y.x;
    C.y.y = A.y.x*B.x.y + A.y.y*B.y.y;

    return C;
}




