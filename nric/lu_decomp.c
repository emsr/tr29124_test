
#include <math.h>



#include "nric.h"



/******************************************************************************************************************

                lu_decomp

        This routine is a double precision version of the routine ludcmp() on p. 43 of NRiC

        Given an n*n matrix a[1..n][1..n], this routine replaces it by the LU (Lower-triangular Upper-triangular) 
    decomposition od a rowwise permutation of itself.  a[][] and n are input.  a[][] is output, index[] is an 
    output vector which row permutation effected by the partial pivoting; d is output as the parity of the row 
    permutation

******************************************************************************************************************/


void  lu_decomp( double **a, int n, int *index, double *parity )

{

    int  i, imax, j, k;
    double  big, dummy, sum, temp;
    double  *vv;

    const double TINY = 1.0e-20;

    vv = dvector( 1, n);
    *parity = 1.0;

    /*
     *    Loop over rows to get the implicit scaling information.    
     */
    for ( i = 1; i <= n; ++i ) {

        big = 0.0;
        for ( j = 1; j <= n; ++j ) {

            temp = fabs( a[i][j]);
            if ( temp > big ) big = temp;
        }
        if ( big == 0.0 ) nrerror( "Singular matrix in routine lu_decomp.");

        /*
         *    Save the scaling.
         */
        vv[i] = 1.0/big;
    }

    /*
     *    This is the loop over columns of Crout's method.
     */
    for ( j = 1; j <= n; ++j ) {

        for ( i = 1; i < j; ++i ) {

            sum = a[i][j];
            for ( k = 1; k < i; ++k ) sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
        }

        /*
         *     Initialize for the search for the largest pivot point.
         */
        big = 0.0;
        for ( i = j; i <= n; ++i ) {

            sum = a[i][j];
            for ( k = 1; k < j; ++k ) sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
            dummy = vv[i]*fabs( sum );
            if ( dummy >= big ) {

                big = dummy;
                imax = i;
            }
        }

        /*
         *     Interchang rows if required.
         */
        if ( j != imax ) {

            for ( k = 1; k <= n; ++k ) {

                dummy = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dummy;
            }

            /*
             *    Change parity.
             */
            *parity = -*parity;

            /*
             *    Interchange the scale factor.
             */
            vv[imax] = vv[j];
        }
        index[j] = imax;
        if ( a[j][j] == 0.0 ) a[j][j] = TINY;

        /*
         *     Now finally devide by the pivot element
         */
        if ( j != n) {

            dummy = 1.0/( a[j][j]);
            for ( i = j + 1; i <= n; ++i ) a[i][j] *= dummy;
        }
    }    /*    Go back for the next column in the reduction.    */
    free_dvector( vv, 1, n);
}




/**************************************************************************************************************************

                LU_BACKSUBSTITUTION

        This routine is a double precision version of the routine on p. 44 of NRiC.

        This routine solves the set of n linear equations a.x=b.  Here a[1..n][1..n] is input, not as the matrix a but as 
    its LU decomposition, determined by the routine lu_decomp().  b[1..n] is input as the right hand side vector b 
    and returns with the left-hand solution vector x.  a, n, and index are not modified by this routine and can be left 
    in place for successive calls with different right hand sides b[1..n].  This routine takes into account the 
    possibility that b will begin with a lot of zeros so that it is efficient for use in matrix inversion.

**************************************************************************************************************************/


void lu_backsub( double **a, int n, int *index, double *b)

{

    int  i, ii=0, ip, j;
    double  sum;

    /*
     *    When ii is set to a posative value, it will become the index of the first nonvanishing element of b[1..n].
     *    We now do the forward substitution.  The only new wrinkle is to unsramble the permutation as we go.
     */
    for ( i = 1; i <= n; ++i ) {

        ip = index[i];
        sum = b[ip];
        b[ip] = b[i];
        if ( ii ) for ( j = ii; j <= i - 1; ++j ) sum -= a[i][j]*b[j];
        else if ( sum ) ii = i;
        b[i] = sum;
    }


    /*
     *    Now do the backsubstitution.
     */

    for ( i = n; i >= 1; i--) {

        sum = b[i];
        for ( j = i + 1; j <= n; ++j ) sum -= a[i][j]*b[j];
        b[i] = sum/a[i][i];
    }
}



/*
 *    Improves a solution vector x of the linear set A.x = b.  The matrix a and the
 *    LU decomposition of a alud (with its row permutation vector index) and the
 *    right-hand side vector are input along with the solution vector x.  
 *    The solution vector x is improved and modified on output.
 */
void lu_improve( double **a, double **alud, int n, int *index, double *b, double *x ) {

    int i, j;
    double *r;

    r = dvector( 1, n );

    for ( i = 1; i <= n; ++i ) {
        r[i] = -b[i];
        for ( j = 1; j <= n; ++j ) r[i] += a[i][j]*x[j];
    }
    lu_backsub( alud, n, index, r );
    for ( i = 1; i <= n; ++i ) x[i] -= r[i];

    free_dvector( r, 1, n );
}



void lu_invert( double **alud, int n, int *index, double **ainv ) {

    int i, j;
    double *col;

    col = dvector( 1, n );

    for ( j = 1; j <= n; ++j ) {
        for ( i = 1; i <= n; ++i ) col[i] = 0.0;
        col[j] = 1.0;
        lu_backsub( alud, n, index, col );
        for ( i = 1; i <= n; ++i ) ainv[i][j] = col[i];
    }

    free_dvector( col, 1, n );
}


/*
 *    Compute determinant of LU decomposed matrix.
 */
double lu_determinant( double **alud, int n, double parity ) {

    int i;
    double det = parity;

    for ( i = 1; i <= n; ++i ) det *= alud[i][i];
    return det;
}


/*
 *    Compute trace of LU decomposed matrix.
 */
double lu_trace( double **alud, int n ) {

    int i, j;
    double trace = 0.0;

    for ( i = 1; i <= n; ++i ) {
        trace += alud[i][i];
        for ( j = i-1; j >= 1; --j ) trace += alud[i][j]*alud[j][i];
    }
    return trace;
}

