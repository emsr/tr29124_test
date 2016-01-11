
#include <math.h>
#include <stdlib.h>

#include "nric.h"


/************************************************************************************
 ************************************************************************************
 **
 **    MIXED routines.
 **
 ************************************************************************************
 ************************************************************************************/

dmatrix32 cart_vec3vec2( dvector3 U, dvector2 V ) {
    dmatrix32 m;

    m.x.x = U.x*V.x;
    m.x.y = U.x*V.y;
    m.y.x = U.y*V.x;
    m.y.y = U.y*V.y;
    m.z.x = U.z*V.x;
    m.z.y = U.z*V.y;

    return m;
}


dmatrix23 cart_vec2vec3( dvector2 U, dvector3 V ) {
    dmatrix23 m;

    m.x.x = U.x*V.x;
    m.x.y = U.x*V.y;
    m.x.z = U.x*V.z;
    m.y.x = U.y*V.x;
    m.y.y = U.y*V.y;
    m.y.z = U.y*V.z;

    return m;
}

dvector3 mat32vec2_dv( dmatrix32 A, dvector2 V ) {
    dvector3 v;

    v.x = A.x.x*V.x + A.x.y*V.y;
    v.y = A.y.x*V.x + A.y.y*V.y;
    v.z = A.z.x*V.x + A.z.y*V.y;

    return v;
}

dvector2 mat23vec3_dv( dmatrix23 A, dvector3 V ) {
    dvector2 v;

    v.x = A.x.x*V.x + A.x.y*V.y + A.x.z*V.z;
    v.y = A.y.x*V.x + A.y.y*V.y + A.y.z*V.z;

    return v;
}

dmatrix32 trans_mat23_dv( dmatrix23 A ) {
    dmatrix32 C;

    C.x.x = A.x.x;
    C.x.y = A.y.x;
    C.y.x = A.x.y;
    C.y.y = A.y.y;
    C.z.x = A.x.z;
    C.z.y = A.y.z;

    return C;
}

dmatrix23 trans_mat32_dv( dmatrix32 A ) {
    dmatrix23 C;

    C.x.x = A.x.x;
    C.x.y = A.y.x;
    C.x.z = A.z.x;
    C.y.x = A.x.y;
    C.y.y = A.y.y;
    C.y.z = A.z.y;

    return C;
}


dmatrix3 mat32mat23_dv( dmatrix32 A, dmatrix23 B ) {
    dmatrix3 C;

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


dmatrix32 mat33mat32_dv( dmatrix3 A, dmatrix32 B ) {
    dmatrix32 C;

    C.x.x = A.x.x*B.x.x + A.x.y*B.y.x + A.x.z*B.z.x;
    C.x.y = A.x.x*B.x.y + A.x.y*B.y.y + A.x.z*B.z.y;
    C.y.x = A.y.x*B.x.x + A.y.y*B.y.x + A.y.z*B.z.x;
    C.y.y = A.y.x*B.x.y + A.y.y*B.y.y + A.y.z*B.z.y;
    C.z.x = A.z.x*B.x.x + A.z.y*B.y.x + A.z.z*B.z.x;
    C.z.y = A.z.x*B.x.y + A.z.y*B.y.y + A.z.z*B.z.y;

    return C;
}


dmatrix32 mat32mat22_dv( dmatrix32 A, dmatrix2 B ) {
    dmatrix32 C;

    C.x.x = A.x.x*B.x.x + A.x.y*B.y.x;
    C.x.y = A.x.x*B.x.y + A.x.y*B.y.y;
    C.y.x = A.y.x*B.x.x + A.y.y*B.y.x;
    C.y.y = A.y.x*B.x.y + A.y.y*B.y.y;
    C.z.x = A.z.x*B.x.x + A.z.y*B.y.x;
    C.z.y = A.z.x*B.x.y + A.z.y*B.y.y;

    return C;
}

dmatrix2 mat23mat32_dv( dmatrix23 A, dmatrix32 B ) {
    dmatrix2 C;

    C.x.x = A.x.x*B.x.x + A.x.y*B.y.x + A.x.z*B.z.x;
    C.x.y = A.x.x*B.x.y + A.x.y*B.y.y + A.x.z*B.z.y;
    C.y.x = A.y.x*B.x.x + A.y.y*B.y.x + A.y.z*B.z.x;
    C.y.y = A.y.x*B.x.y + A.y.y*B.y.y + A.y.z*B.z.y;

    return C;
}


dmatrix23 mat22mat23_dv( dmatrix2 A, dmatrix23 B ) {
    dmatrix23 C;

    C.x.x = A.x.x*B.x.x + A.x.y*B.y.x;
    C.x.y = A.x.x*B.x.y + A.x.y*B.y.y;
    C.x.z = A.x.x*B.x.z + A.x.y*B.y.z;
    C.y.x = A.y.x*B.x.x + A.y.y*B.y.x;
    C.y.y = A.y.x*B.x.y + A.y.y*B.y.y;
    C.y.z = A.y.x*B.x.z + A.y.y*B.y.z;

    return C;
}


dmatrix23 mat23mat33_dv( dmatrix23 A, dmatrix3 B ) {
    dmatrix23 C;

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


dvector3 *dv3vector( int nl, int nh ) {

    dvector3 *v;

    v = ( dvector3 * ) malloc( ( unsigned long ) ( nh - nl + 1 )*sizeof( dvector3 ) );
    if ( !v ) nrerror( "Allocation failure in dv3vector." );
    v -= nl;

    return  v;
}

void  free_dv3vector( dvector3 *v, int nl, int nh ) {
    free( ( void * ) ( v + nl ) );
}


dvector2 *dv2vector( int nl, int nh ) {

    dvector2 *v;

    v = ( dvector2 * ) malloc( ( unsigned long ) ( nh - nl + 1 )*sizeof( dvector2 ) );
    if ( !v ) nrerror( "Allocation failure in dv2vector." );
    v -= nl;

    return  v;
}

void  free_dv2vector( dvector2 *v, int nl, int nh ) {
    free( ( void * ) ( v + nl ) );
}


double *convert_dvector_dv3( dvector3 a ) {

    double *v;
    v = dvector( 1, 3 );
    v[1] = a.x;
    v[2] = a.y;
    v[3] = a.z;
    return v;
}

void free_convert_dvector_dv3( double *v ) {

    free_dvector( v, 1, 3 );
}


double **convert_dmatrix_dv3( dmatrix3 a ) {

    double **m;
    m = dmatrix( 1, 3, 1, 3 );
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

void free_convert_dmatrix_dv3( double **m ) {

    free_dmatrix( m, 1, 3, 1, 3 );
}


dvector3 mk_dvector_dv3( double *v ) {

    dvector3 a;
    a.x = v[1];
    a.y = v[2];
    a.z = v[3];
    return a;
}


dmatrix3 mk_dmatrix_dv3( double **m ) {

    dmatrix3 a;
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


dvector3 mk_dv3( double x, double y, double z ) {

    dvector3 c;
    c.x = x;
    c.y = y;
    c.z = z;
    return c;
}



dvector3 add_dv3( dvector3 a, dvector3 b ) {

    dvector3 c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    return c;
}




dvector3 sub_dv3( dvector3 a, dvector3 b ) {

    dvector3 c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    return c;
}




double mod_dv3( dvector3 a ) {

    double max, myx, rx, ry, rz;

    max = ((rz = fabs(a.z)) >= (myx = ((ry = fabs(a.y)) >= (rx = fabs(a.x)) ? ry : rx)) ? rz : myx );

    if ( max == 0.0 ) return 0.0;
    else {
        rx /= max;
        ry /= max;
        rz /= max;
        return max*sqrt( rx*rx + ry*ry + rz*rz );
    }
}





dvector3 norm_dv3( dvector3 a ) {

    double max, myx, rx, ry, rz;

    max = ((rz = fabs(a.z)) >= (myx = ((ry = fabs(a.y)) >= (rx = fabs(a.x)) ? ry : rx)) ? rz : myx );

    if ( max == 0.0 ) {
        a.x = a.y = a.z = 0.0;
    } else {
        rx /= max;
        ry /= max;
        rz /= max;
        max *= sqrt( rx*rx + ry*ry + rz*rz );

        a.x /= max;
        a.y /= max;
        a.z /= max;
    }
    return a;
}




double dot_dv3( dvector3 a, dvector3 b ) {

    double c;
    c = a.x*b.x + a.y*b.y + a.z*b.z;
    return c;
}




double az_dv3( dvector3 a ) {

    double c;
    c = atan2( a.y, a.x );
    return c;
}




double el_dv3( dvector3 a ) {

    double c;
    c = atan2( a.z, sqrt( a.x*a.x + a.y*a.y ) );
    return c;
}




dvector3 cross_dv3( dvector3 a, dvector3 b ) {

    dvector3 c;
    c.x = a.y*b.z - a.z*b.y;
    c.y = a.z*b.x - a.x*b.z;
    c.z = a.x*b.y - a.y*b.x;
    return c;
}




dvector3 rmul_dv3( double a, dvector3 b ) {

    dvector3 c;
    c.x = a*b.x;
    c.y = a*b.y;
    c.z = a*b.z;
    return c;
}




dmatrix3 cart_dv3( dvector3 a, dvector3 b ) {

    dmatrix3 c;
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



double trace_dv3( dmatrix3 a ) { return a.x.x + a.y.y + a.z.z; }

double det_dv3( dmatrix3 A ) {

     return A.x.x*(A.y.y*A.z.z - A.y.z*A.z.y)
          + A.x.y*(A.y.z*A.z.x - A.y.x*A.z.z)
          + A.x.z*(A.y.x*A.z.y - A.y.y*A.z.x);
}

dmatrix3 inv_dv3( dmatrix3 A ) {

    double det;
    dmatrix3 a;

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


dmatrix3 trans_dv3( dmatrix3 A ) {

    dmatrix3 a;

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



/************************************************************************************
 ************************************************************************************
 **
 **    VECTOR2 routines.
 **
 ************************************************************************************
 ************************************************************************************/



double *convert_dvector_dv2( dvector2 a ) {

    double *v;
    v = dvector( 1, 2 );
    v[1] = a.x;
    v[2] = a.y;
    return v;
}


void free_convert_dvector_dv2( double *v ) {

    free_dvector( v, 1, 2 );
}


double **convert_dmatrix_dv2( dmatrix2 a ) {

    double **m;
    m = dmatrix( 1, 2, 1, 2 );
    m[1][1] = a.x.x;
    m[1][2] = a.x.y;
    m[2][1] = a.y.x;
    m[2][2] = a.y.y;
    return m;
}

void free_convert_dmatrix_dv2( double **m ) {

    free_dmatrix( m, 1, 2, 1, 2 );
}


dvector2 mk_dvector_dv2( double *v ) {

    dvector2 a;
    a.x = v[1];
    a.y = v[2];
    return a;
}


dmatrix2 mk_dmatrix_dv2( double **m ) {

    dmatrix2 a;
    a.x.x = m[1][1];
    a.x.y = m[1][2];
    a.y.x = m[2][1];
    a.y.y = m[2][2];
    return a;
}


dvector2 mk_dv2( double x, double y ) {

    dvector2 c;
    c.x = x;
    c.y = y;
    return c;
}



dvector2 add_dv2( dvector2 a, dvector2 b ) {

    dvector2 c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    return c;
}




dvector2 sub_dv2( dvector2 a, dvector2 b ) {

    dvector2 c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    return c;
}




double mod_dv2( dvector2 a ) {

    double max, rx, ry;

    max = ((ry = fabs(a.y)) >= (rx = fabs(a.x)) ? ry : rx);

    if ( max == 0.0 ) return 0.0;
    else {
        rx /= max;
        ry /= max;
        return max*sqrt( rx*rx + ry*ry );
    }
}




dvector2 norm_dv2( dvector2 a ) {

    double max, rx, ry;

    max = ((ry = fabs(a.y)) >= (rx = fabs(a.x)) ? ry : rx);

    if ( max == 0.0 ) {
        a.x = a.y = 0.0;
    } else {
        rx /= max;
        ry /= max;
        max *= sqrt( rx*rx + ry*ry );

        a.x /= max;
        a.y /= max;
    }
    return a;
}




double dot_dv2( dvector2 a, dvector2 b ) {

    double c;
    c = a.x*b.x + a.y*b.y ;
    return c;
}




double az_dv2( dvector2 a ) {

    double c;
    c = atan2( a.y, a.x );
    return c;
}



dvector2 rmul_dv2( double a, dvector2 b ) {

    dvector2 c;
    c.x = a*b.x;
    c.y = a*b.y;
    return c;
}




dmatrix2 cart_dv2( dvector2 a, dvector2 b ) {

    dmatrix2 c;
    c.x.x = a.x*b.x;
    c.x.y = a.x*b.y;
    c.y.x = a.y*b.x;
    c.y.y = a.y*b.y;
    return c;
}


dmatrix2 trans_dv2( dmatrix2 A ) {

    dmatrix2 a;

    a.x.x = A.x.x;
    a.x.y = A.y.x;
    a.y.x = A.x.y;
    a.y.y = A.y.y;

    return a;
}


double trace_dv2( dmatrix2 A ) { return A.x.x + A.y.y; }


double det_dv2( dmatrix2 A ) { return A.x.x*A.y.y - A.x.y*A.y.x; }


dmatrix2 inv_dv2( dmatrix2 A ) {

    double det;
    dmatrix2 a;

    det = A.x.x*A.y.y - A.x.y*A.y.x;

    a.x.x = A.y.y/det;
    a.y.x = -A.y.x/det;
    a.x.y = A.x.y/det;
    a.y.y = -A.x.x/det;

    return a;
}



dvector3 matvec_dv3( dmatrix3 A, dvector3 X ) {

    dvector3 V;

    V.x = A.x.x*X.x + A.x.y*X.y + A.x.z*X.z;
    V.y = A.y.x*X.x + A.y.y*X.y + A.y.z*X.z;
    V.z = A.z.x*X.x + A.z.y*X.y + A.z.z*X.z;

    return V;
}



dvector3 vecmat_dv3( dvector3 Y, dmatrix3 A ) {

    dvector3 V;

    V.x = Y.x*A.x.x + Y.y*A.y.x + Y.z*A.z.x;
    V.y = Y.x*A.x.y + Y.y*A.y.y + Y.z*A.z.y;
    V.z = Y.x*A.x.z + Y.y*A.y.z + Y.z*A.z.z;

    return V;
}



double vecmatvec_dv3( dvector3 Y, dmatrix3 A, dvector3 X ) {

    return (Y.x*A.x.x + Y.y*A.y.x + Y.z*A.z.x)*X.x
         + (Y.x*A.x.y + Y.y*A.y.y + Y.z*A.z.y)*X.y
         + (Y.x*A.x.z + Y.y*A.y.z + Y.z*A.z.z)*X.z;
}



dmatrix3 matmat_dv3( dmatrix3 A, dmatrix3 B ) {

    dmatrix3 C;

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







dmatrix2 matmat_dv2( dmatrix2 A, dmatrix2 B ) {

    dmatrix2 C;

    C.x.x = A.x.x*B.x.x + A.x.y*B.y.x;
    C.x.y = A.x.x*B.x.y + A.x.y*B.y.y;
    C.y.x = A.y.x*B.x.x + A.y.y*B.y.x;
    C.y.y = A.y.x*B.x.y + A.y.y*B.y.y;

    return C;
}





dvector2 matvec_dv2( dmatrix2 A, dvector2 X ) {

    dvector2 V;

    V.x = A.x.x*X.x + A.x.y*X.y;
    V.y = A.y.x*X.x + A.y.y*X.y;

    return V;
}



dvector2 vecmat_dv2( dvector2 Y, dmatrix2 A ) {

    dvector2 V;

    V.x = Y.x*A.x.x + Y.y*A.y.x;
    V.y = Y.x*A.x.y + Y.y*A.y.y;

    return V;
}



double vecmatvec_dv2( dvector2 Y, dmatrix2 A, dvector2 X ) {

    return (Y.x*A.x.x + Y.y*A.y.x)*X.x
         + (Y.x*A.x.y + Y.y*A.y.y)*X.y;
}






