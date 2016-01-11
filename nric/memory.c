#include <stdio.h>
#include <malloc.h>

#include "nric.h"



/*****************************************************************************

                UCVECTOR

    Allocates an unsigned char vector with range [nl..nh].

*****************************************************************************/

unsigned char *ucvector( int nl, int nh ) {

    unsigned char *v;

    v = ( unsigned char * ) malloc( ( unsigned int ) ( nh - nl + 1 )*sizeof( unsigned char ) );
    if ( !v ) nrerror( "Allocation failure in ucvector." );
    v -= nl;

    return v;
}




/*****************************************************************************

                FREE_UCVECTOR

    Frees an unsigned char vector allocated with ucvector().

*****************************************************************************/

void  free_ucvector( unsigned char *v, int nl, int nh ) {

    free( ( void * ) ( v + nl ) );
}




/*****************************************************************************

                CVECTOR

    Allocates an char vector with range [nl..nh].

*****************************************************************************/

char *cvector( int nl, int nh ) {

    char *v;

    v = ( char * ) malloc( ( unsigned int ) ( nh - nl + 1 )*sizeof( char ) );
    if ( !v ) nrerror( "Allocation failure in cvector." );
    v -= nl;

    return v;
}




/*****************************************************************************

                FREE_CVECTOR

    Frees an char vector allocated with cvector().

*****************************************************************************/

void  free_cvector( char *v, int nl, int nh ) {

    free( ( void * ) ( v + nl ) );
}



/*****************************************************************************

                CMATRIX

    Allocates a char matrix with range [nrl..nrh][ncl..nch] and returns
  a pointer to an array of pointers to rows.

*****************************************************************************/

char **cmatrix( int nrl, int nrh, int ncl, int nch ) {

    int i;
    char **m;

    /*
     *    Allocate pointers to rows.
     */
    m = ( char ** ) malloc( ( unsigned int ) ( nrh - nrl + 1 )*sizeof( char * ) );
    if ( !m ) nrerror( "Allocation failure 1 in cmatrix." );
    m -= nrl;


    /*
     *    Allocate rows and set pointers to them.
     */
    for ( i = nrl; i <= nrh; ++i) {

        m[i] = ( char * ) malloc( ( unsigned int ) ( nch - ncl + 1 )*sizeof( char ) );
        if ( !m[i] ) nrerror( "Allocation failure 2 in cmatrix." );
        m[i] -= ncl;
    }

    return  m;
}



/*****************************************************************************

                FREE_CMATRIX

    Frees a char matrix m, allocated with cmatrix.

*****************************************************************************/

void  free_cmatrix( char **m, int nrl, int nrh, int ncl, int nch ) {

    int i;

    for ( i = nrh; i >= nrl; --i ) free( ( void * ) ( m[i] + ncl ) );
    free( ( void * ) ( m + nrl ) );
}




/*****************************************************************************

                UIVECTOR

    Allocates an unsigned integer vector with range [nl..nh].

*****************************************************************************/

unsigned int *uivector( int nl, int nh ) {

    unsigned int *v;

    v = ( unsigned int * ) malloc( ( unsigned int ) ( nh - nl + 1 )*sizeof( unsigned int ) );
    if ( !v ) nrerror( "Allocation failure in uivector." );
    v -= nl;

    return v;
}




/*****************************************************************************

                FREE_UIVECTOR

    Frees an integer vector allocated with uivector().

*****************************************************************************/

void  free_uivector( unsigned int *v, int nl, int nh ) {

    free( ( void * ) ( v + nl ) );
}





/*****************************************************************************

                IVECTOR

    Allocates an integer vector with range [nl..nh].

*****************************************************************************/

int *ivector( int nl, int nh ) {

    int *v;

    v = ( int * ) malloc( ( unsigned int ) ( nh - nl + 1 )*sizeof( int ) );
    if ( !v ) nrerror( "Allocation failure in ivector." );
    v -= nl;

    return v;
}




/*****************************************************************************

                FREE_IVECTOR

    Frees an integer vector allocated with ivector().

*****************************************************************************/

void  free_ivector( int *v, int nl, int nh ) {

    free( ( void * ) ( v + nl ) );
}




/*****************************************************************************

                LVECTOR

    Allocates a long vector with range [nl..nh].

*****************************************************************************/

long *lvector( int nl, int nh ) {

    long *v;

    v = ( long * ) malloc( ( unsigned int ) ( nh - nl + 1 )*sizeof( long ) );
    if ( !v ) nrerror( "Allocation failure in lvector." );
    v -= nl;

    return v;
}




/*****************************************************************************

                FREE_LVECTOR

    Frees a long vector allocated with lvector().

*****************************************************************************/

void  free_lvector( long *v, int nl, int nh ) {

    free( ( void * ) ( v + nl ) );
}




/*****************************************************************************

                ULVECTOR

    Allocates a long vector with range [nl..nh].

*****************************************************************************/

unsigned long *ulvector( int nl, int nh ) {

    unsigned long *v;

    v = ( unsigned long * ) malloc( ( unsigned int ) ( nh - nl + 1 )*sizeof( unsigned long ) );
    if ( !v ) nrerror( "Allocation failure in ulvector." );
    v -= nl;

    return v;
}




/*****************************************************************************

                FREE_ULVECTOR

    Frees a long vector allocated with ulvector().

*****************************************************************************/

void  free_ulvector( unsigned long *v, int nl, int nh ) {

    free( ( void * ) ( v + nl ) );
}





/*****************************************************************************

                IMATRIX

    Allocates a int matrix with range [nrl..nrh][ncl..nch] and returns
  a pointer to an array of pointers to rows.

*****************************************************************************/

int **imatrix( int nrl, int nrh, int ncl, int nch ) {

    int i;
    int **m;

    /*
     *    Allocate pointers to rows.
     */
    m = ( int ** ) malloc( ( unsigned int ) ( nrh - nrl + 1 )*sizeof( int * ) );
    if ( !m ) nrerror( "Allocation failure 1 in imatrix." );
    m -= nrl;


    /*
     *    Allocate rows and set pointers to them.
     */
    for ( i = nrl; i <= nrh; ++i ) {
        m[i] = ( int * ) malloc( ( unsigned int ) ( nch - ncl + 1 )*sizeof( int ) );
        if ( !m[i] ) nrerror( "Allocation failure 2 in imatrix." );
        m[i] -= ncl;
    }

    return m;
}



/*****************************************************************************

                FREE_IMATRIX

    Frees a int matrix m, allocated with imatrix.

*****************************************************************************/

void  free_imatrix( int **m, int nrl, int nrh, int ncl, int nch ) {

    int i;

    for ( i = nrh; i >= nrl; --i ) free( ( void * ) ( m[i] + ncl ) );
    free( ( void * ) ( m + nrl ) );
}




/*****************************************************************************

                FVECTOR

    Allocates a float vector with range [nl..nh].

*****************************************************************************/

float *fvector( int nl, int nh ) {

    float *v;

    v = ( float * ) malloc( ( unsigned int ) ( nh - nl + 1 )*sizeof( float ) );
    if ( !v ) nrerror( "Allocation failure in fvector." );
    v -= nl;

    return  v;
}



/*****************************************************************************

                FREE_FVECTOR

    Frees a float vector allocated with fvector().

*****************************************************************************/

void  free_fvector( float *v, int nl, int nh ) {

    free( ( void * ) ( v + nl ) );
}



/*****************************************************************************

                FMATRIX

    Allocates a float matrix with range [nrl..nrh][ncl..nch] and returns
  a pointer to an array of pointers to rows.

*****************************************************************************/

float **fmatrix( int nrl, int nrh, int ncl, int nch ) {

    int i;
    float **m;

    /*
     *    Allocate pointers to rows.
     */
    m = ( float ** ) malloc( ( unsigned int ) ( nrh - nrl + 1 )*sizeof( float * ) );
    if ( !m ) nrerror( "Allocation failure 1 in fmatrix." );
    m -= nrl;


    /*
     *    Allocate rows and set pointers to them.
     */
    for ( i = nrl; i <= nrh; ++i) {

        m[i] = ( float * ) malloc( ( unsigned int ) ( nch - ncl + 1 )*sizeof( float ) );
        if ( !m[i] ) nrerror( "Allocation failure 2 in fmatrix." );
        m[i] -= ncl;
    }

    return  m;
}



/*****************************************************************************

                FREE_FMATRIX

    Frees a float matrix m, allocated with fmatrix.

*****************************************************************************/

void  free_fmatrix( float **m, int nrl, int nrh, int ncl, int nch ) {

    int i;

    for ( i = nrh; i >= nrl; --i ) free( ( void * ) ( m[i] + ncl ) );
    free( ( void * ) ( m + nrl ) );
}




/*****************************************************************************

                CONVERT_FMATRIX

    Alocates a float matrix m, that points to the float matrix a, declared 
  in the standard C manner as a[nrow][ncol] where nrow = nrh - nrl + 1 and
  ncol = nch - ncl + 1.

    The routine should be called with the address &a[0][0] as the first 
  argument.

*****************************************************************/

float **convert_fmatrix( float *a, int nrl, int nrh, int ncl, int nch ) {

    int i, j, nrow, ncol;
    float **m;


    nrow = nrh - nrl + 1;
    ncol = nch - ncl + 1;

    m = ( float ** ) malloc( ( unsigned int ) nrow*sizeof( float * ) );
    if ( !m ) nrerror( "Allocation failure 1 in convert_fmatrix." );
    m -= nrl;

    for ( i = 1, j = nrl; i <= nrow - 1; ++i, ++j)  {

        m[j] = a + ncol*i - ncl;
        if ( !m[i] ) nrerror( "Allocation failure 2 in convert_fmatrix." );
        m[i] -= ncl;
    }

    return  m;
}




/*****************************************************************

                FREE_CONVERT_FMATRIX

    Free a matrix allocated with convert_dmatrix().

*****************************************************************/

void  free_convert_fmatrix( float **m, int nrl, int nrh, int ncl, int nch ) {

    free( ( void * ) ( m + nrl ) );
}





/*****************************************************************************

                DVECTOR

    Allocates a double vector with range [nl..nh].

*****************************************************************************/

double *dvector( int nl, int nh ) {

    double *v;

    v = ( double * ) malloc( ( unsigned int ) ( nh - nl + 1 )*sizeof( double ) );
    if ( !v ) nrerror( "Allocation failure in dvector." );
    v -= nl;

    return  v;
}



/*****************************************************************************

                FREE_DVECTOR

    Frees a double vector allocated with dvector().

*****************************************************************************/

void  free_dvector( double *v, int nl, int nh ) {
    free( ( void * ) ( v + nl ) );
}



/*****************************************************************************

                DMATRIX

    Allocates a double matrix with range [nrl..nrh][ncl..nch] and returns
  a pointer to an array of pointers to rows.

*****************************************************************************/

double **dmatrix( int nrl, int nrh, int ncl, int nch ) {

    int i;
    double **m;

    /*
     *    Allocate pointers to rows.
     */
    m = ( double ** ) malloc( ( unsigned int ) ( nrh - nrl + 1 )*sizeof( double * ) );
    if ( !m ) nrerror( "Allocation failure 1 in dmatrix." );
    m -= nrl;


    /*
     *    Allocate rows and set pointers to them.
     */
    for ( i = nrl; i <= nrh; ++i) {

        m[i] = ( double * ) malloc( ( unsigned int ) ( nch - ncl + 1 )*sizeof( double ) );
        if ( !m[i] ) nrerror( "Allocation failure 2 in dmatrix." );
        m[i] -= ncl;
    }

    return  m;
}



/*****************************************************************************

                FREE_DMATRIX

    Frees a double matrix m, allocated with dmatrix.

*****************************************************************************/

void  free_dmatrix( double **m, int nrl, int nrh, int ncl, int nch ) {

    int i;

    for ( i = nrh; i >= nrl; --i ) free( ( void * ) ( m[i] + ncl ) );
    free( ( void * ) ( m + nrl ) );
}




/*****************************************************************************

                CONVERT_DMATRIX

    Alocates a double matrix m, that points to the double matrix a, declared 
  in the standard C manner as a[nrow][ncol] where nrow = nrh - nrl + 1 and
  ncol = nch - ncl + 1.

    The routine should be called with the address &a[0][0] as the first 
  argument.

*****************************************************************************/

double **convert_dmatrix( double *a, int nrl, int nrh, int ncl, int nch ) {

    int i, j, nrow, ncol;
    double **m;

    nrow = nrh - nrl + 1;
    ncol = nch - ncl + 1;

    m = ( double ** ) malloc( ( unsigned int ) nrow*sizeof( double * ) );
    if ( !m ) nrerror( "Allocation failure 1 in convert_dmatrix." );
    m -= nrl;

    for ( i = 1, j = nrl; i <= nrow - 1; ++i, ++j )  {

        m[j] = a + ncol*i - ncl;
        if ( !m[i] ) nrerror( "Allocation failure 2 in convert_dmatrix." );
        m[i] -= ncl;
    }

    return  m;
}




/*****************************************************************

                FREE_CONVERT_DMATRIX

    Free a matrix allocated with convert_dmatrix().

*****************************************************************/

void  free_convert_dmatrix( double **m, int nrl, int nrh, int ncl, int nch ) {
    free( ( void * ) ( m + nrl ) );
}


