
#include <math.h>


#include "nric.h"


  /*
   *    Linear equation solution by Gauss-Jordan elimination. a[1..n][1..n] is an input matrix.
   *    b[1..n][1..m] is the input matrix of m right-hand side vectors.
   *    On output, a is replaced by its inverse, and b is replaced by the corresponding set of
   *    solution vectors.
   */
  void
  gauss_jordan(double** a, int n, double** b, int m)
  {
    int icol, irow;
    double dum;

    int* indxc = ivector(1, n);
    int* indxr = ivector(1, n);
    int* ipiv = ivector(1, n);

    for (int j = 1; j <= n; ++j)
      ipiv[j] = 0;
    for (int i = 1; i <= n; ++i)
      {
	/*
	 * Loop over columns to be reduced.
	 */
	double big = 0.0;
	for (int j = 1; j <= n; ++j)
	  {
	    /*
	     * Loop over rows looking for the pivot elements.
	     */
	    if (ipiv[j] != 1)
	      for (int k = 1; k <= n; ++k)
		{
		  if (ipiv[k] == 0)
		    {
		      if (fabs(a[j][k]) >= big)
			{
			  big = fabs(a[j][k]);
			  irow = j;
			  icol = k;
			}
		    }
		  else if (ipiv[k] > 1)
		      nrerror("Singular matrix in gauss_jordan.");
		}
	  }
	++(ipiv[icol]);
printf("ipiv[icol]=%d\n", ipiv[icol]);

	/*
	 * With the pivot elements in hand, we swap rows to put the pivot elements on the diagonal.
	 * The columns are not physically moved, only relabeled: ipiv[i], the column of the
	 * original ith pivot element, is the ith column that is reduced, while indxr[i] is the
	 * row in which that pivot element was originally located.
	 * If indxr[i] != indxc[i] an implied column interchange.
	 */
	if (irow != icol)
	  {
	    for (int l = 1; l <= n; ++l)
	      dswap(&a[irow][l], &a[icol][l]);
	    for (int l = 1; l <= m; ++l)
	      dswap(&b[irow][l], &b[icol][l]);
	  }
	indxr[i] = irow;
	indxc[i] = icol;
printf("indxr[i]=%d\n", indxr[i]);
printf("indxc[i]=%d\n", indxc[i]);
	if (a[icol][icol] == 0.0)
	  nrerror("Singular matrix error 2 in gauss_jordan.");
	double pivinv = 1.0 / a[icol][icol];
	a[icol][icol] = 1.0;
	for (int l = 1; l <= n; ++l)
	  a[icol][l] *= pivinv;
	for (int l = 1; l <= m; ++l)
	  b[icol][l] *= pivinv;
	for (int ll = 1; ll <= n; ++ll)
	  {
	    if (ll != icol)
	      {
		dum = a[ll][icol];
		a[ll][icol] = 0.0;
		for (int l = 1; l <= n; ++l)
		  a[ll][l] -= a[icol][l] * dum;
		for (int l = 1; l <= m; ++l)
		  b[ll][l] -= b[icol][l] * dum;
	      }
	  }
      }
    for (int l = n; l >= 1; --l)
      {
	if (indxr[l] != indxc[l])
	  for (int k = 1; k <= n; ++k)
	    dswap(&a[k][indxr[l]], &a[k][indxc[l]]);
      }

    free_ivector(ipiv, 1, n);
    free_ivector(indxr, 1, n);
    free_ivector(indxc, 1, n);
  }



