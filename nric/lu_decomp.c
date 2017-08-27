
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


  void
  lu_decomp(double** a, int n, int* index, double* parity)
  {
    const double TINY = 1.0e-20;

    double* vv = dvector(1, n);
    *parity = 1.0;

    /* Loop over rows to get the implicit scaling information. */
    for (int i = 1; i <= n; ++i)
      {
	double big = 0.0;
	for (int j = 1; j <= n; ++j)
	  {
	    double temp = fabs(a[i][j]);
	    if (temp > big)
	      big = temp;
	  }
	if (big == 0.0)
	    nrerror("Singular matrix in routine lu_decomp.");

	/*
	 * Save the scaling.
	 */
	vv[i] = 1.0 / big;
printf("big=%f\n", big);
printf("scale[i]=%f\n", vv[i]);
      }

    /* This is the loop over columns of Crout's method. */
    for (int j = 1; j <= n; ++j)
      {
	for (int i = 1; i < j; ++i)
	  {
	    double sum = a[i][j];
	    for (int k = 1; k < i; ++k)
	      sum -= a[i][k] * a[k][j];
	    a[i][j] = sum;
printf("a[i][j]=%f\n", a[i][j]);
	  }

	/* Initialize for the search for the largest pivot point. */
	int imax = 0;
	double big = 0.0;
	for (int i = j; i <= n; ++i)
	  {
	    double sum = a[i][j];
	    for (int k = 1; k < j; ++k)
	      sum -= a[i][k] * a[k][j];
	    a[i][j] = sum;
printf("a[i][j]=%f\n", a[i][j]);
	    double dummy = vv[i] * fabs(sum);
	    if (dummy >= big)
	      {
		big = dummy;
		imax = i;
	      }
	  }
printf("big=%f\n", big);
printf("imax=%d\n", imax);

	/* Interchang rows if required. */
	if (j != imax)
	  {
	    for (int k = 1; k <= n; ++k)
	      {
		double dummy = a[imax][k];
		a[imax][k] = a[j][k];
		a[j][k] = dummy;
	      }

	    /* Change parity. */
	    *parity = -*parity;

	    /* Interchange the scale factor. */
	    vv[imax] = vv[j];
	  }
	index[j] = imax;
printf("index[j]=%d\n", index[j]);
	if (a[j][j] == 0.0)
	    a[j][j] = TINY;

	/* Now finally divide by the pivot element. */
	if (j != n)
	  {
	    double dummy = 1.0 / a[j][j];
	    for (int i = j + 1; i <= n; ++i)
{
	      a[i][j] *= dummy;
printf("a[i][j]=%f\n", a[i][j]);
}
	  }
      } /* Go back for the next column in the reduction. */

    free_dvector(vv, 1, n);
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


  void
  lu_backsub(double** a, int n, int* index, double* b)
  {
    int  ii=0;

    /*
     * When ii is set to a posative value, it will become the index of the first nonvanishing element of b[1..n].
     * We now do the forward substitution.  The only new wrinkle is to unsramble the permutation as we go.
     */
    for (int i = 1; i <= n; ++i)
      {
	int ip = index[i];
	double sum = b[ip];
	b[ip] = b[i];
	if (ii)
	  for (int j = ii; j <= i - 1; ++j)
	    sum -= a[i][j] * b[j];
	else if (sum)
	  ii = i;
	b[i] = sum;
      }


    /* Now do the backsubstitution. */
    for (int i = n; i >= 1; i--)
      {
	double sum = b[i];
	for (int j = i + 1; j <= n; ++j)
	  sum -= a[i][j] * b[j];
	b[i] = sum/a[i][i];
      }
  }



/*
 *    Improves a solution vector x of the linear set A.x = b.  The matrix a and the
 *    LU decomposition of a alud (with its row permutation vector index) and the
 *    right-hand side vector are input along with the solution vector x.
 *    The solution vector x is improved and modified on output.
 */
  void
  lu_improve(double** a, double** alud, int n, int* index, double* b, double* x)
  {
    double* r = dvector(1, n);

    for (int i = 1; i <= n; ++i)
      {
	r[i] = -b[i];
	for (int j = 1; j <= n; ++j)
	  r[i] += a[i][j]*x[j];
      }
    lu_backsub(alud, n, index, r);
    for (int i = 1; i <= n; ++i)
      x[i] -= r[i];

    free_dvector(r, 1, n);
  }



  void
  lu_invert(double** alud, int n, int* index, double** ainv)
  {
    double* col = dvector(1, n);

    for (int j = 1; j <= n; ++j)
      {
	for (int i = 1; i <= n; ++i)
	  col[i] = 0.0;
	col[j] = 1.0;
	lu_backsub(alud, n, index, col);
	for (int i = 1; i <= n; ++i)
	  ainv[i][j] = col[i];
      }

    free_dvector(col, 1, n);
  }


  /*
   * Compute determinant of LU decomposed matrix.
   */
  double
  lu_determinant(double** alud, int n, double parity)
  {
    double det = parity;
    for (int i = 1; i <= n; ++i)
      det *= alud[i][i];
    return det;
  }


  /*
   * Compute trace of LU decomposed matrix.
   */
  double
  lu_trace(double** alud, int n)
  {
    double trace = 0.0;
    for (int i = 1; i <= n; ++i)
      {
	trace += alud[i][i];
	for (int j = i-1; j >= 1; --j)
	  trace += alud[i][j]*alud[j][i];
      }
    return trace;
  }

