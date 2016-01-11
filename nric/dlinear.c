
#include <math.h>

#include "nric.h"



/**************************************************************************************
 *
 *    This is a library of utility routines for basic linear algebra of vectors
 *    and matrices.  Routines are available for copying vectors and matrices, for
 *    multiplying vectors, matrix and vector, vector and matrix, and matrix and matrix.
 *    Routines are available for adding, subtracting, vectors and matrices and for
 *    multiplying vectors and matrices by scalars.
 *
 *    All the routines here require that the matrices and vectors be allocated
 *    with malloc or using the dvector dmatrix allocation routines.
 *
 *************************************************************************************/
/*****************************************************************************************
 *
 *    Double versions.
 *
 *****************************************************************************************/



void dcopyvec( double *a, int nrowlo, int nrowhi, double *b ) {
    int i;
    for ( i = nrowlo; i <= nrowhi; ++i ) b[i] = a[i];
}



void dcopymat( double **a, int nrowlo, int nrowhi, int ncollo, int ncolhi, double **b ) {
    int i, j;
    for ( i = nrowlo; i <= nrowhi; ++i ) for ( j = ncollo; j <= ncolhi; ++j ) b[i][j] = a[i][j];
}



/*
 *    Input is a[nlo..nhi].
 *    Output is c.
 */
void dmodvec( double *a, int nlo, int nhi, double *c ) {
    int k;
    double max, *r;

    r = dvector( nlo, nhi );

    max = r[nlo] = fabs(a[nlo]);
    for ( k = nlo+1; k <= nhi; ++k ) if ( ( r[k] = fabs(a[k]) ) > max ) max = r[k];

    if ( max == 0.0 ) *c = 0.0;
    else {
        for ( *c = 0.0, k = nlo; k <= nhi; ++k ) {
            r[k] /= max;
            *c += r[k]*r[k];
        }
        *c = max*sqrt(*c);
    }

    free_dvector( r, nlo, nhi );
}



/*
 *    Input is a[nrowlo..nrowhi] and b[ncollo..ncolhi].
 *    Output is c[nrowlo..nrowhi][ncollo..ncolhi].
 */
void dcartvecvec( double *a, double *b, int nrowlo, int nrowhi, int ncollo, int ncolhi, double **c ) {
    int i, j;
    for ( i = nrowlo; i <= nrowhi; ++i ) for ( j = ncollo; j <= ncolhi; ++j ) c[i][j] = a[i]*b[j];
}



/*
 *    Input is a[nlo..nhi] and b[nlo..nhi].
 *    Output is c.
 */
void dvecvec( double *a, double *b, int nlo, int nhi, double *c ) {
    int k;
    for ( *c = 0.0, k = nlo; k <= nhi; ++k ) *c += a[k]*b[k];
}



/*
 *    Input is a[nrowlo..nrowhi][nlo..nhi], and b[nlo..nhi].
 *    Output is c[nrowlo..nrowhi].
 */
void dmatvec( double **a, double *b, int nrowlo, int nrowhi, int nlo, int nhi, double *c ) {
    int i, k;
    for ( i = nrowlo; i <= nrowhi; ++i ) {
        c[i] = 0.0;
        for ( k = nlo; k <= nhi; ++k ) c[i] += a[i][k]*b[k];
    }
}


/*
 *    Input is b[nlo..nhi], a[nlo..nhi][ncollo..ncolhi]
 *    Output is c[ncollo..ncolhi].
 */
void dvecmat( double *b, double **a, int nlo, int nhi, int ncollo, int ncolhi, double *c ) {
    int k, j;
    for ( j = ncollo; j <= ncolhi; ++j ) {
        c[j] = 0.0;
        for ( k = nlo; k <= nhi; ++k ) c[j] += b[k]*a[k][j];
    }
}


/*
 *    Input is v[n1lo..n1hi], a[n1lo..n1hi][n2lo..n2hi], and u[n2lo..n2hi].
 *    Output is c.
 */
void dvecmatvec( double *v, double **a, double *u, int n1lo, int n1hi, int n2lo, int n2hi, double *c ) {
    int i, j;
    for ( *c = 0.0, i = n1lo; i <= n1hi; ++i ) for ( j = n2lo; j <= n2hi; ++j ) *c += v[i]*a[i][j]*u[j];
}




/*
 *    Input is a[nrowlo..nrowhi][nlo..nhi], b[nlo..nhi][ncollo..ncolhi]
 *    Output is c[nrowlo][nrowhi][ncollo..ncolhi].
 */
void dmatmat( double **a, double **b, int nrowlo, int nrowhi, int nlo, int nhi, int ncollo, int ncolhi, double **c ) {
    int i, j, k;
    for ( i = nrowlo; i <= nrowhi; ++i ) {
        for ( j = ncollo; j <= ncolhi; ++j ) {
            c[i][j] = 0.0;
            for ( k = nlo; k <= nhi; ++k ) { c[i][j] += a[i][k]*b[k][j]; }
        }
    }
}


/*
 *    Input is a[nrowlo..nrowhi], b[nrowlo..nrowhi]
 *    Output is c[nrowlo..nrowhi].
 */
void daddvec( double *a, double *b, int nrowlo, int nrowhi, double *c ) {
    int i;
    for ( i = nrowlo; i <= nrowhi; ++i ) c[i] = a[i] + b[i];
}


/*
 *    Input is a[nrowlo..nrowhi], b[nrowlo..nrowhi]
 *    Output is c[nrowlo][nrowhi].
 */
void dsubvec( double *a, double *b, int nrowlo, int nrowhi, double *c ) {
    int i;
    for ( i = nrowlo; i <= nrowhi; ++i ) c[i] = a[i] - b[i];
}


/*
 *    Input is a, b[nrowlo..nrowhi]
 *    Output is c[nrowlo][nrowhi].
 */
void dmulvec( double a, double *b, int nrowlo, int nrowhi, double *c ) {
    int i;
    for ( i = nrowlo; i <= nrowhi; ++i ) c[i] = a*b[i];
}



/*
 *    Input is a[nrowlo..nrowhi][ncollo..ncolhi], b[nrowlo..nrowhi][ncollo..ncolhi]
 *    Output is c[nrowlo][nrowhi][ncollo..ncolhi].
 */
void daddmat( double **a, double **b, int nrowlo, int nrowhi, int ncollo, int ncolhi, double **c ) {
    int i, j;
    for ( i = nrowlo; i <= nrowhi; ++i ) for ( j = ncollo; j <= ncolhi; ++j ) c[i][j] = a[i][j] + b[i][j];
}


/*
 *    Input is a[nrowlo..nrowhi][ncollo..ncolhi], b[nrowlo..nrowhi][ncollo..ncolhi]
 *    Output is c[nrowlo][nrowhi][ncollo..ncolhi].
 */
void dsubmat( double **a, double **b, int nrowlo, int nrowhi, int ncollo, int ncolhi, double **c ) {
    int i, j;
    for ( i = nrowlo; i <= nrowhi; ++i ) for ( j = ncollo; j <= ncolhi; ++j ) c[i][j] = a[i][j] - b[i][j];
}


/*
 *    Input is a, b[nrowlo..nrowhi][ncollo..ncolhi]
 *    Output is c[nrowlo][nrowhi][ncollo..ncolhi].
 */
void dmulmat( double a, double **b, int nrowlo, int nrowhi, int ncollo, int ncolhi, double **c ) {
    int i, j;
    for ( i = nrowlo; i <= nrowhi; ++i ) for ( j = ncollo; j <= ncolhi; ++j ) c[i][j] = a*b[i][j];
}



/*
 *    Input is a square matrix a[nrowlo..nrowhi][ncollo..ncolhi],
 *    Output is the transpose atrans[ncollo..ncolhi][nrowlo..nrowhi].
 */
void dtransmat( double **a, int nrowlo, int nrowhi, int ncollo, int ncolhi, double **atrans ) {
    int i, j;
    for ( i = ncollo; i <= ncolhi; ++i ) for ( j = nrowlo; j <= nrowhi; ++j ) atrans[i][j] = a[j][i];
}



/*
 *    Input is a[nlo..nhi][nlo..nhi],
 *    Output is the symmetric part of a asymm[nlo..nhi][nlo..nhi],
 *    and the antisymmetric part of a aanti[nlo..nhi][nlo..nhi].
 */
void dsymmat( double **a, int nlo, int nhi, double **asymm, double **aanti ) {
    int i, j;
    for ( i = nlo; i <= nhi; ++i ) {
        for ( j = i; j <= nhi; ++j ) {
            asymm[j][i] =  ( asymm[i][j] = ( a[i][j] + a[j][i] )/2.0 );
            aanti[j][i] = -( aanti[i][j] = ( a[i][j] - a[j][i] )/2.0 );
        }
    }
}



/*
 *    Input is a[nlo..nhi][nlo..nhi], a square matrix.
 *    Output is the trace of the matrix.
 */
void dtracemat( double **a, int nlo, int nhi, double *trace ) {
    int i;
    *trace = 0.0;
    for ( i = nlo; i <= nhi; ++i ) *trace += a[i][i];
}



/*
 *    Gramm-Schmidt orthogonalization.
 */
void dgramm_schmidt( double **a, int nlo, int nhi, double **aortho ) {

    int i, j;
    double dot, norm, *ortho;

    ortho = dvector( nlo, nhi );

    /*
     *    Normalize first column vector as the first basis vector.
     *    Initialize the orthoganal vector.
     */
    for ( norm = 0.0, i = nlo; i <= nhi; ++i ) norm += a[i][nlo];
    for ( norm = sqrt(norm), i = nlo; i <= nhi; ++i ) ortho[i] = aortho[i][nlo] /= norm;

    /*
     *    Loop through the remaining columns.
     *    Subtract components of previous column vectors from each new vector
     *    and normalize.
     */
    for ( j = nlo+1; j <= nhi; ++j ) {
        for ( dot = 0.0, i = nlo; i <= nhi; ++i ) dot += ortho[i]*a[i][j];
        for ( norm = 0.0, i = nlo; i <= nhi; ++i ) norm += ( aortho[i][j] = a[i][j] - dot*ortho[i] );
        for ( norm = sqrt(norm), i = nlo; i <= nhi; ++i ) ortho[i] += ( aortho[i][j] /= norm );
    }

    free_dvector( ortho, nlo, nhi );
}



