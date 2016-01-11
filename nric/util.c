

#include <math.h>


#include "nric.h"


/*
 *    Returns the sign of the argument.
 */

double dsgn( double a ) { return a >= 0.0 ? 1.0 : -1.0; }

float fsgn( float a ) { return a >= 0.0 ? 1.0 : -1.0; }

int isgn( int a ) { return a >= 0 ? 1 : -1; }



/*
 *    Returns the first argument with the sign of the second argument.
 */

double dsign( double a, double b ) {
    return a == 0.0 ? 0.0 : ( b >= 0.0 ? fabs(a) : -fabs(a) );
}

float fsign( float a, float b ) {
    return a == 0.0 ? 0.0 : ( b >= 0.0 ? fabsf(a) : -fabsf(a) );
}

int isign( int a, int b ) {
    return a == 0.0 ? 0 : ( b >= 0 ? (a >= 0 ? a : -a) : -(a >= 0 ? a : -a) );
}


/*
 *    Returns the minimum argument.
 */

double dmin( double a, double b ) {
    return a < b ? a : b;
}

//float fmin( float a, float b ) {
//    return a < b ? a : b;
//}

int imin( int a, int b ) {
    return a < b ? a : b;
}


/*
 *    Returns the minimum argument.
 */

double dmin3( double a, double b, double c ) {
    double bc;
    return a < (bc=(b < c ? b : c)) ? a : bc;
}

float fmin3( float a, float b, float c ) {
    float bc;
    return a < (bc=(b < c ? b : c)) ? a : bc;
}

int imin3( int a, int b, int c ) {
    int bc;
    return a < (bc=(b < c ? b : c)) ? a : bc;
}



/*
 *    Returns the maximum argument.
 */

double dmax( double a, double b ) {
    return b >= a ? b : a;
}

//float fmax( float a, float b ) {
//    return b >= a ? b : a;
//}

int imax( int a, int b ) {
    return b >= a ? b : a;
}


/*
 *    Returns the maximum of three arguments.
 */

double dmax3( double a, double b, double c ) {
    double ba;
    return c >= (ba=(b >= a ? b : a)) ? c : ba;
}

float fmax3( float a, float b, float c ) {
    float ba;
    return c >= (ba=(b >= a ? b : a)) ? c : ba;
}

int imax3( int a, int b, int c ) {
    int ba;
    return c >= (ba=(b >= a ? b : a)) ? c : ba;
}





/*
 *    Swaps the value of two numbers.
 */

void dswap( double *a, double *b ) {
    double temp = *a;
    *a = *b;
    *b = temp;
}

void fswap( float *a, float *b ) {
    float temp = *a;
    *a = *b;
    *b = temp;
}

void iswap( int *a, int *b ) {
    int temp = *a;
    *a = *b;
    *b = temp;
}



/*
 *    Permutes the values of three numbers in a right-handed sense.
 */

void drpermute3( double *a, double *b, double *c ) {
    double temp = *c;
    *c = *b;
    *b = *a;
    *a = temp;
}

void frpermute3( float *a, float *b, float *c ) {
    float temp = *c;
    *c = *b;
    *b = *a;
    *a = temp;
}

void irpermute3( int *a, int *b, int *c ) {
    int temp = *c;
    *c = *b;
    *b = *a;
    *a = temp;
}



/*
 *    Permutes the values of three numbers in a left-handed sense.
 */

void dlpermute3( double *a, double *b, double *c ) {
    double temp = *a;
    *a = *b;
    *b = *c;
    *c = temp;
}

void flpermute3( float *a, float *b, float *c ) {
    float temp = *a;
    *a = *b;
    *b = *c;
    *c = temp;
}

void ilpermute3( int *a, int *b, int *c ) {
    int temp = *a;
    *a = *b;
    *b = *c;
    *c = temp;
}



/*
 *    Shifts the values of three numbers in a right-handed sense.
 *    The old right-most element is overwritten.
 *    The left-most element is set to zero.
 */

void drshift3( double *a, double *b, double *c ) {
    *c = *b;
    *b = *a;
    *a = 0.0;
}

void frshift3( float *a, float *b, float *c ) {
    *c = *b;
    *b = *a;
    *a = 0.0;
}

void irshift3( int *a, int *b, int *c ) {
    *c = *b;
    *b = *a;
    *a = 0;
}



/*
 *    Shifts the values of three numbers in a left-handed sense.
 *    The old left-most element is overwritten.
 *    The right-most element is set to zero.
 */

void dlshift3( double *a, double *b, double *c ) {
    *a = *b;
    *b = *c;
    *c = 0.0;
}

void flshift3( float *a, float *b, float *c ) {
    *a = *b;
    *b = *c;
    *c = 0.0;
}

void ilshift3( int *a, int *b, int *c ) {
    *a = *b;
    *b = *c;
    *c = 0;
}



/*
 *    Shifts the values of four numbers in a right-handed sense.
 *    The old right-most element is overwritten.
 *    The left-most element is set to zero.
 */

void drshift4( double *a, double *b, double *c, double *d) {
    *d = *c;
    *c = *b;
    *b = *a;
    *a = 0.0;
}

void frshift4( float *a, float *b, float *c, float *d) {
    *d = *c;
    *c = *b;
    *b = *a;
    *a = 0.0;
}

void irshift4( int *a, int *b, int *c, int *d) {
    *d = *c;
    *c = *b;
    *b = *a;
    *a = 0;
}



/*
 *    Shifts the values of three numbers in a left-handed sense.
 *    The old left-most element is overwritten.
 *    The right-most element is set to zero.
 */

void dlshift4( double *a, double *b, double *c, double *d) {
    *a = *b;
    *b = *c;
    *c = *d;
    *d = 0.0;
}

void flshift4( float *a, float *b, float *c, float *d) {
    *a = *b;
    *b = *c;
    *c = *d;
    *d = 0.0;
}

void ilshift4( int *a, int *b, int *c, int *d) {
    *a = *b;
    *b = *c;
    *c = *d;
    *d = 0;
}



/*
 *    These return the square of thier arguments.
 */

double dsqr( double a ) { return a == 0.0 ? 0.0 : a*a; }

float fsqr( float a ) { return a == 0.0 ? 0.0 : a*a; }

int isqr( int a ) { return a == 0 ? 0 : a*a; }



/*
 *    These return the cube of thier arguments.
 */

double dcub( double a ) { return a == 0.0 ? 0.0 : a*a*a; }

float fcub( float a ) { return a == 0.0 ? 0.0 : a*a*a; }

int icub( int a ) { return a == 0 ? 0 : a*a*a; }



/*
 *    These return the 4th power of thier arguments.
 */

double dpow4( double a ) { return a == 0.0 ? 0.0 : (a *= a *= a); }

float fpow4( float a ) { return a == 0.0 ? 0.0 : (a *= a *= a); }

int ipow4( int a ) { return a == 0 ? 0 : (a *= a *= a); }



/*
 *    These return the 5th power of thier arguments.
 */

double dpow5( double a ) { double A = a; return a == 0.0 ? 0.0 : A*(a *= a *= a); }

float fpow5( float a ) { float A = a; return a == 0.0 ? 0.0 : A*(a *= a *= a); }

int ipow5( int a ) { int A = a; return a == 0 ? 0 : A*(a *= a *= a); }


/*
 *    Computes sqrt(a^2 + b^2) without destructive overflow or underflow.
 */

double dpythag( double a, double b ) {
    double at, bt, ct;
    return ((at=fabs(a)) > (bt=fabs(b)) ? (ct=bt/at, at*sqrt(1.0 + ct*ct)) : (bt ? (ct=at/bt, bt*sqrt(1.0 + ct*ct)) : 0.0));
}

float fpythag( float a, float b ) {
    float at, bt, ct;
    return ((at=fabsf(a)) > (bt=fabsf(b)) ? (ct=bt/at, at*sqrt(1.0 + ct*ct)) : (bt ? (ct=at/bt, bt*sqrt(1.0 + ct*ct)) : 0.0));
}

